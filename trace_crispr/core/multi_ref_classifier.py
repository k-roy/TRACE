"""
Alignment-based classification for CRISPR editing outcomes.

This module classifies reads based on alignment to a multi-reference FASTA
containing WT + all HDR template variants. Classification is based on:
1. Which reference the read aligns to (WT vs HDR variants)
2. Alignment quality (mismatches, soft clips)
3. Presence of indels near the cut site

Author: Kevin R. Roy
Date: 2026-02-19
"""

import logging
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pysam

logger = logging.getLogger(__name__)


@dataclass
class Indel:
    """Represents an insertion or deletion."""
    position: int  # 0-based position in reference
    length: int    # Positive = deletion, negative = insertion
    sequence: str  # Inserted/deleted sequence if available

    @property
    def is_insertion(self) -> bool:
        return self.length < 0

    @property
    def is_deletion(self) -> bool:
        return self.length > 0

    @property
    def size(self) -> int:
        return abs(self.length)


@dataclass
class AlignmentClassification:
    """Classification result for a single read."""
    outcome: str  # 'WT', 'HDR_PERFECT', 'HDR_IMPERFECT', 'NHEJ', 'LARGE_DEL', 'UNMAPPED'
    aligned_ref: str  # Reference name (WT, HDR_bc1, etc.)
    template_id: Optional[str] = None  # Barcode ID if HDR
    mismatches: int = 0
    indels_near_cut: List[Indel] = field(default_factory=list)
    best_aligner: str = ""
    alignment_score: int = 0
    soft_clip_5: int = 0  # Soft clip at 5' end
    soft_clip_3: int = 0  # Soft clip at 3' end
    read_count: int = 1  # For collapsed reads with count in header


@dataclass
class ClassificationSummary:
    """Summary of classification results for a sample."""
    total_reads: int = 0
    wt_count: int = 0
    hdr_perfect_count: int = 0
    hdr_imperfect_count: int = 0
    nhej_count: int = 0
    large_del_count: int = 0
    unmapped_count: int = 0

    # Per-barcode HDR counts
    hdr_by_barcode: Dict[str, int] = field(default_factory=dict)

    # Expected barcode tracking
    expected_barcode: Optional[str] = None
    expected_barcode_count: int = 0

    @property
    def hdr_total_count(self) -> int:
        return self.hdr_perfect_count + self.hdr_imperfect_count

    @property
    def wt_rate(self) -> float:
        return self.wt_count / self.total_reads if self.total_reads > 0 else 0.0

    @property
    def hdr_rate(self) -> float:
        return self.hdr_total_count / self.total_reads if self.total_reads > 0 else 0.0

    @property
    def nhej_rate(self) -> float:
        return self.nhej_count / self.total_reads if self.total_reads > 0 else 0.0

    @property
    def large_del_rate(self) -> float:
        return self.large_del_count / self.total_reads if self.total_reads > 0 else 0.0

    @property
    def expected_barcode_rate(self) -> float:
        return self.expected_barcode_count / self.total_reads if self.total_reads > 0 else 0.0


class MultiRefClassifier:
    """
    Classify reads based on alignment to multi-reference FASTA.

    Classification logic:
    - Aligned to HDR + 0 mismatches, no indels → HDR_PERFECT
    - Aligned to HDR + 1-5 mismatches, no indels → HDR_IMPERFECT
    - Aligned to HDR + >5 mismatches OR indels → Inspect further (may be NHEJ)
    - Aligned to WT + indels within cut_site_window → NHEJ
    - Aligned to WT + deletion ≥50bp spanning cut → LARGE_DEL
    - Aligned to WT + no indels, ≤5 mismatches → WT
    - Unmapped → UNMAPPED
    """

    def __init__(
        self,
        cut_sites: Dict[str, int],
        cut_site_window: int = 10,
        max_mismatches_imperfect: int = 5,
        large_del_threshold: int = 50,
    ):
        """
        Initialize the classifier.

        Args:
            cut_sites: Dict mapping reference name to cut site position (0-based)
            cut_site_window: bp around cut site for NHEJ detection
            max_mismatches_imperfect: Max mismatches for imperfect HDR
            large_del_threshold: Minimum deletion size for LARGE_DEL classification
        """
        self.cut_sites = cut_sites
        self.cut_site_window = cut_site_window
        self.max_mismatches_imperfect = max_mismatches_imperfect
        self.large_del_threshold = large_del_threshold

    def classify_from_bams(
        self,
        bam_paths: Dict[str, Path],
        expected_barcode: Optional[str] = None,
    ) -> Tuple[Dict[str, AlignmentClassification], ClassificationSummary]:
        """
        Parse BAM files from multiple aligners and classify reads.

        Args:
            bam_paths: Dict mapping aligner name to BAM path
            expected_barcode: Expected HDR barcode for this sample (optional)

        Returns:
            Tuple of (per-read classifications, summary statistics)
        """
        # Collect alignments from all BAMs
        all_alignments = defaultdict(list)  # read_name -> [(aligner, alignment)]

        for aligner, bam_path in bam_paths.items():
            if not bam_path.exists():
                logger.warning(f"BAM file not found: {bam_path}")
                continue

            try:
                with pysam.AlignmentFile(str(bam_path), "rb") as bam:
                    for read in bam:
                        all_alignments[read.query_name].append((aligner, read))
            except Exception as e:
                logger.error(f"Error reading {bam_path}: {e}")
                continue

        # Classify each read
        classifications = {}
        summary = ClassificationSummary(expected_barcode=expected_barcode)

        for read_name, alignments in all_alignments.items():
            # Get read count from header (for collapsed reads: seq_N_count=X)
            read_count = self._extract_read_count(read_name)

            # Select best alignment
            best_aligner, best_alignment = self._select_best_alignment(alignments)

            # Classify
            classification = self._classify_alignment(
                best_alignment, best_aligner, read_count
            )
            classifications[read_name] = classification

            # Update summary
            self._update_summary(summary, classification, expected_barcode)

        return classifications, summary

    def _extract_read_count(self, read_name: str) -> int:
        """Extract read count from collapsed read header."""
        # Format: seq_N_count=X or just seq_N
        if "_count=" in read_name:
            try:
                return int(read_name.split("_count=")[1].split("_")[0])
            except (ValueError, IndexError):
                pass
        return 1

    def _select_best_alignment(
        self,
        alignments: List[Tuple[str, pysam.AlignedSegment]]
    ) -> Tuple[str, pysam.AlignedSegment]:
        """
        Select the best alignment from multiple aligners.

        Selection criteria (in order):
        1. Mapped alignments preferred over unmapped
        2. Higher alignment score (AS tag)
        3. Fewer mismatches (NM tag)
        """
        # Filter to mapped alignments
        mapped = [(a, r) for a, r in alignments if not r.is_unmapped]

        if not mapped:
            # Return any alignment (all unmapped)
            return alignments[0]

        # Sort by alignment score (higher is better), then by edit distance (lower is better)
        def score_key(item):
            aligner, read = item
            as_score = read.get_tag('AS') if read.has_tag('AS') else 0
            nm = read.get_tag('NM') if read.has_tag('NM') else 100
            return (-as_score, nm)

        mapped.sort(key=score_key)
        return mapped[0]

    def _classify_alignment(
        self,
        alignment: pysam.AlignedSegment,
        aligner: str,
        read_count: int,
    ) -> AlignmentClassification:
        """Classify a single alignment."""

        if alignment.is_unmapped:
            return AlignmentClassification(
                outcome='UNMAPPED',
                aligned_ref='',
                best_aligner=aligner,
                read_count=read_count,
            )

        ref_name = alignment.reference_name
        is_hdr = ref_name.startswith('HDR_')
        template_id = ref_name.replace('HDR_', '') if is_hdr else None

        # Get alignment metrics
        mismatches = alignment.get_tag('NM') if alignment.has_tag('NM') else 0
        as_score = alignment.get_tag('AS') if alignment.has_tag('AS') else 0

        # Parse CIGAR for indels and soft clips
        indels, soft_clip_5, soft_clip_3 = self._parse_cigar(alignment)

        # Get cut site for this reference
        cut_site = self.cut_sites.get(ref_name, None)

        # Find indels near cut site
        indels_near_cut = []
        if cut_site is not None:
            cut_window_start = cut_site - self.cut_site_window
            cut_window_end = cut_site + self.cut_site_window

            for indel in indels:
                # Check if indel overlaps with cut site window
                indel_end = indel.position + indel.size if indel.is_deletion else indel.position
                if (indel.position <= cut_window_end and indel_end >= cut_window_start):
                    indels_near_cut.append(indel)

        # Classify based on reference and alignment quality
        outcome = self._determine_outcome(
            is_hdr=is_hdr,
            mismatches=mismatches,
            indels_near_cut=indels_near_cut,
            all_indels=indels,
            has_any_indels=len(indels) > 0,
        )

        return AlignmentClassification(
            outcome=outcome,
            aligned_ref=ref_name,
            template_id=template_id,
            mismatches=mismatches,
            indels_near_cut=indels_near_cut,
            best_aligner=aligner,
            alignment_score=as_score,
            soft_clip_5=soft_clip_5,
            soft_clip_3=soft_clip_3,
            read_count=read_count,
        )

    def _parse_cigar(
        self,
        alignment: pysam.AlignedSegment
    ) -> Tuple[List[Indel], int, int]:
        """
        Parse CIGAR string for indels and soft clips.

        Returns:
            Tuple of (list of indels, 5' soft clip, 3' soft clip)
        """
        indels = []
        soft_clip_5 = 0
        soft_clip_3 = 0

        if alignment.cigartuples is None:
            return indels, soft_clip_5, soft_clip_3

        ref_pos = alignment.reference_start
        query_pos = 0

        for i, (op, length) in enumerate(alignment.cigartuples):
            # CIGAR operations:
            # 0 = M (match/mismatch)
            # 1 = I (insertion)
            # 2 = D (deletion)
            # 4 = S (soft clip)
            # 5 = H (hard clip)

            if op == 4:  # Soft clip
                if i == 0:
                    soft_clip_5 = length
                else:
                    soft_clip_3 = length
                query_pos += length
            elif op == 1:  # Insertion
                seq = alignment.query_sequence[query_pos:query_pos+length] if alignment.query_sequence else ""
                indels.append(Indel(position=ref_pos, length=-length, sequence=seq))
                query_pos += length
            elif op == 2:  # Deletion
                indels.append(Indel(position=ref_pos, length=length, sequence=""))
                ref_pos += length
            elif op == 0:  # Match/mismatch
                ref_pos += length
                query_pos += length
            elif op == 5:  # Hard clip - doesn't consume query
                pass

        return indels, soft_clip_5, soft_clip_3

    def _determine_outcome(
        self,
        is_hdr: bool,
        mismatches: int,
        indels_near_cut: List[Indel],
        all_indels: List[Indel],
        has_any_indels: bool = False,
    ) -> str:
        """Determine the classification outcome."""

        has_indels_near_cut = len(indels_near_cut) > 0

        # Check for large deletions (applies to both WT and HDR alignments)
        for indel in all_indels:
            if indel.is_deletion and indel.size >= self.large_del_threshold:
                return 'LARGE_DEL'

        if is_hdr:
            # HDR classification - aligned to an HDR template
            # HDR templates have the guide/PAM modified, so we don't check cut site
            # Just check alignment quality (mismatches and indels)
            if has_any_indels:
                # Any indels in HDR alignment suggests imperfect repair
                if mismatches <= self.max_mismatches_imperfect:
                    return 'HDR_IMPERFECT'
                else:
                    # Too many mismatches + indels - ambiguous
                    return 'NHEJ'
            elif mismatches == 0:
                return 'HDR_PERFECT'
            elif mismatches <= self.max_mismatches_imperfect:
                return 'HDR_IMPERFECT'
            else:
                # Too many mismatches - likely poor alignment or NHEJ
                return 'NHEJ'
        else:
            # WT classification - aligned to wild-type reference
            if has_indels_near_cut:
                return 'NHEJ'
            elif mismatches <= self.max_mismatches_imperfect:
                return 'WT'
            else:
                # Too many mismatches to be clean WT
                return 'NHEJ'

    def _update_summary(
        self,
        summary: ClassificationSummary,
        classification: AlignmentClassification,
        expected_barcode: Optional[str],
    ):
        """Update summary statistics with a classification."""
        count = classification.read_count
        summary.total_reads += count

        if classification.outcome == 'WT':
            summary.wt_count += count
        elif classification.outcome == 'HDR_PERFECT':
            summary.hdr_perfect_count += count
            if classification.template_id:
                summary.hdr_by_barcode[classification.template_id] = \
                    summary.hdr_by_barcode.get(classification.template_id, 0) + count
                if expected_barcode and classification.template_id == expected_barcode:
                    summary.expected_barcode_count += count
        elif classification.outcome == 'HDR_IMPERFECT':
            summary.hdr_imperfect_count += count
            if classification.template_id:
                summary.hdr_by_barcode[classification.template_id] = \
                    summary.hdr_by_barcode.get(classification.template_id, 0) + count
                if expected_barcode and classification.template_id == expected_barcode:
                    summary.expected_barcode_count += count
        elif classification.outcome == 'NHEJ':
            summary.nhej_count += count
        elif classification.outcome == 'LARGE_DEL':
            summary.large_del_count += count
        elif classification.outcome == 'UNMAPPED':
            summary.unmapped_count += count


def classify_sample(
    bam_paths: Dict[str, Path],
    cut_sites: Dict[str, int],
    expected_barcode: Optional[str] = None,
    cut_site_window: int = 10,
) -> ClassificationSummary:
    """
    Convenience function to classify a single sample.

    Args:
        bam_paths: Dict mapping aligner name to BAM path
        cut_sites: Dict mapping reference name to cut site position
        expected_barcode: Expected HDR barcode for this sample
        cut_site_window: bp around cut site for NHEJ detection

    Returns:
        ClassificationSummary with all statistics
    """
    classifier = MultiRefClassifier(
        cut_sites=cut_sites,
        cut_site_window=cut_site_window,
    )

    _, summary = classifier.classify_from_bams(
        bam_paths=bam_paths,
        expected_barcode=expected_barcode,
    )

    return summary
