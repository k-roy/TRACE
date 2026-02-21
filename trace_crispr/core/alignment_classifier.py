"""
Triple aligner-based classification for NHEJ and large deletion detection.

This module wraps BWA-MEM, BBMap, and minimap2 to provide comprehensive
CRISPR editing outcome classification including:
- HDR detection (alignment-based validation)
- NHEJ insertions and deletions
- Large deletions (≥50bp)

Author: Kevin R. Roy
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
from collections import defaultdict
import logging
import tempfile
import pysam

from ..integrations.aligners import run_triple_alignment, create_reference_fasta
from .scoring import score_alignment, AlignmentScore
from .cigar import get_deletions_from_cigar, get_insertions_from_cigar
from .classification import check_indels_at_cut_site, EditingOutcome

logger = logging.getLogger(__name__)


@dataclass
class AlignmentConfig:
    """Configuration for alignment-based classification."""
    large_deletion_threshold: int = 50  # bp
    cleavage_window_size: int = 12  # bp around cut site
    min_anchor_length: int = 20  # bp for deletion validation
    max_mismatches_20bp: int = 1  # for anchor validation
    max_mismatches_50bp: int = 2  # for anchor validation
    threads: int = 1  # aligner threads (SLURM safety - set to 1 to avoid resource overuse)


@dataclass
class AlignmentClassificationResult:
    """Result from alignment-based classification."""
    outcome: EditingOutcome
    template_id: Optional[str] = None  # For HDR, which template
    is_perfect: bool = False  # For HDR, perfect vs imperfect
    indel_info: Optional[Dict] = None  # For NHEJ, details about indel
    best_aligner: Optional[str] = None  # Which aligner was selected
    alignment_score: Optional[float] = None  # Penalty score from best alignment


class AlignmentClassifier:
    """
    Triple aligner-based classifier for NHEJ and large deletion detection.

    Workflow:
    1. Create FASTQ from sequences
    2. Run BWA-MEM, BBMap, minimap2
    3. Score and select best alignment per sequence
    4. Parse CIGARs for indels
    5. Classify outcomes (HDR, WT, NHEJ, large deletion)
    """

    def __init__(
        self,
        reference: str,
        guide: str,
        hdr_templates: Optional[Dict[str, str]] = None,
        config: Optional[AlignmentConfig] = None,
    ):
        """
        Initialize alignment classifier.

        Args:
            reference: Reference sequence (WT)
            guide: Guide RNA sequence
            hdr_templates: Optional dict of template_id -> HDR sequence
            config: AlignmentConfig object (uses defaults if not provided)
        """
        self.reference = reference
        self.guide = guide
        self.hdr_templates = hdr_templates or {}
        self.config = config or AlignmentConfig()

        # Calculate cut site from guide (3bp upstream of PAM)
        # For SpCas9 with NGG PAM, cut site is ~3bp upstream of PAM (blunt-ended DSB)
        if guide in reference.upper():
            guide_pos = reference.upper().index(guide.upper())
            # PAM is immediately 3' of guide
            pam_pos = guide_pos + len(guide)
            self.cut_site = pam_pos - 3  # Cut site ~3bp UPSTREAM of PAM (CORRECTED)
            logger.info(f"Guide found at position {guide_pos}, PAM at {pam_pos}, cut site at {self.cut_site} (3bp upstream of PAM)")
        else:
            # If guide not in reference, use middle of reference
            self.cut_site = len(reference) // 2
            logger.warning(f"Guide not found in reference, using middle of reference as cut site: {self.cut_site}")

    def classify_sequences(
        self,
        sequences: Set[str],
        temp_dir: Path,
    ) -> Dict[str, AlignmentClassificationResult]:
        """
        Classify a set of unique sequences using triple alignment.

        Args:
            sequences: Set of unique sequences to classify
            temp_dir: Temporary directory for intermediate files

        Returns:
            Dict mapping sequence -> AlignmentClassificationResult
        """
        if not sequences:
            return {}

        logger.info(f"Classifying {len(sequences)} sequences using triple alignment...")

        temp_dir.mkdir(parents=True, exist_ok=True)

        # Step 1: Create reference FASTA
        ref_fasta = temp_dir / "reference.fasta"
        create_reference_fasta(self.reference, ref_fasta, name="reference")
        logger.info(f"Created reference FASTA: {ref_fasta}")

        # Step 2: Create query FASTQ from sequences
        query_fastq = temp_dir / "query.fastq"
        self._write_sequences_to_fastq(sequences, query_fastq)
        logger.info(f"Created query FASTQ with {len(sequences)} sequences: {query_fastq}")

        # Step 3: Run triple alignment
        alignment_dir = temp_dir / "alignments"
        alignment_results = run_triple_alignment(
            r1_fastq=query_fastq,
            r2_fastq=None,  # Single-end
            reference=ref_fasta,
            output_dir=alignment_dir,
            threads=self.config.threads,
        )

        # Step 4: Parse alignments and select best per sequence
        best_alignments = self._select_best_alignments(alignment_results, sequences)
        logger.info(f"Selected best alignments for {len(best_alignments)} sequences")

        # Step 5: Classify based on best alignments
        classifications = self._classify_from_alignments(best_alignments)
        logger.info(f"Classified {len(classifications)} sequences")

        # Log classification summary
        outcome_counts = defaultdict(int)
        for result in classifications.values():
            outcome_counts[result.outcome.value] += 1
        logger.info(f"Alignment classification summary: {dict(outcome_counts)}")

        return classifications

    def _write_sequences_to_fastq(self, sequences: Set[str], output_path: Path):
        """Write sequences to FASTQ with dummy quality scores."""
        with open(output_path, 'w') as f:
            for i, seq in enumerate(sequences, 1):
                # Use sequence index as read ID
                f.write(f"@seq_{i}\n")
                f.write(f"{seq}\n")
                f.write("+\n")
                f.write("I" * len(seq) + "\n")  # Dummy high quality

    def _select_best_alignments(
        self,
        alignment_results: Dict[str, any],
        sequences: Set[str],
    ) -> Dict[str, Tuple[AlignmentScore, pysam.AlignedSegment]]:
        """
        Select best alignment per sequence across all three aligners.

        Args:
            alignment_results: Dict from run_triple_alignment
            sequences: Original sequences (for mapping read IDs)

        Returns:
            Dict mapping sequence -> (AlignmentScore, AlignedSegment)
        """
        # Build sequence lookup by read ID
        seq_list = list(sequences)
        read_id_to_seq = {f"seq_{i+1}": seq for i, seq in enumerate(seq_list)}

        # Parse all alignments and score them
        alignments_by_read = defaultdict(list)

        for aligner_name, aligner_result in alignment_results.items():
            if not aligner_result.success:
                logger.warning(f"{aligner_name} failed, skipping")
                continue

            # Parse BAM
            try:
                with pysam.AlignmentFile(aligner_result.bam_path, 'rb') as bam:
                    for read in bam:
                        # Score this alignment
                        score = score_alignment(read, aligner_name)
                        alignments_by_read[read.query_name].append((score, read))
            except Exception as e:
                logger.error(f"Error parsing BAM from {aligner_name}: {e}")
                continue

        # Select best alignment per read
        best_alignments = {}
        for read_id, alignment_list in alignments_by_read.items():
            if not alignment_list:
                continue

            # Sort by score (lower is better)
            best_score, best_read = min(alignment_list, key=lambda x: x[0])

            # Map back to original sequence
            seq = read_id_to_seq.get(read_id)
            if seq:
                best_alignments[seq] = (best_score, best_read)

        return best_alignments

    def _classify_from_alignments(
        self,
        best_alignments: Dict[str, Tuple[AlignmentScore, pysam.AlignedSegment]],
    ) -> Dict[str, AlignmentClassificationResult]:
        """
        Classify sequences based on their best alignments.

        Classification priority:
        1. NHEJ insertions (insertions in cleavage window)
        2. Large deletions (≥50bp spanning edit site)
        3. NHEJ deletions (deletions in cleavage window)
        4. HDR (if templates provided)
        5. WT (perfect match to reference)
        6. Unclassified
        """
        classifications = {}

        for seq, (score, read) in best_alignments.items():
            if not score.is_aligned:
                classifications[seq] = AlignmentClassificationResult(
                    outcome=EditingOutcome.UNMAPPED
                )
                continue

            # Check for NHEJ indels near cut site
            deletions, insertions = check_indels_at_cut_site(
                read,
                self.cut_site,
                window=self.config.cleavage_window_size,
            )

            # Priority 1: NHEJ insertions
            if insertions:
                classifications[seq] = AlignmentClassificationResult(
                    outcome=EditingOutcome.NHEJ_INSERTION,
                    indel_info={'insertions': insertions},
                    best_aligner=score.aligner,
                    alignment_score=score.penalty_score,
                )
                continue

            # Priority 2: Large deletions
            all_deletions = get_deletions_from_cigar(read, min_size=1)
            for deletion in all_deletions:
                if deletion.size >= self.config.large_deletion_threshold:
                    # Check if deletion spans cut site
                    if deletion.ref_start <= self.cut_site <= deletion.ref_end:
                        classifications[seq] = AlignmentClassificationResult(
                            outcome=EditingOutcome.LARGE_DELETION,
                            indel_info={'deletion': {
                                'start': deletion.ref_start,
                                'end': deletion.ref_end,
                                'size': deletion.size,
                            }},
                            best_aligner=score.aligner,
                            alignment_score=score.penalty_score,
                        )
                        continue

            # Priority 3: NHEJ deletions (small deletions in cleavage window)
            if deletions:
                classifications[seq] = AlignmentClassificationResult(
                    outcome=EditingOutcome.NHEJ_DELETION,
                    indel_info={'deletions': deletions},
                    best_aligner=score.aligner,
                    alignment_score=score.penalty_score,
                )
                continue

            # Priority 4: HDR detection (if templates provided)
            # TODO: Implement HDR signature matching similar to k-mer approach
            # For now, skip HDR detection - k-mer should handle HDR well

            # Priority 5: WT (no indels, good alignment)
            if score.penalty_score == 0:
                classifications[seq] = AlignmentClassificationResult(
                    outcome=EditingOutcome.WILD_TYPE,
                    best_aligner=score.aligner,
                    alignment_score=score.penalty_score,
                )
                continue

            # Default: Unclassified
            classifications[seq] = AlignmentClassificationResult(
                outcome=EditingOutcome.UNCLASSIFIED,
                best_aligner=score.aligner,
                alignment_score=score.penalty_score,
            )

        return classifications
