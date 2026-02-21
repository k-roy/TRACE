"""
Batch processing for multi-template CRISPR editing analysis.

This module provides optimized batch processing for experiments with multiple
HDR templates (e.g., barcoded CRISPR screens). The key optimization is global
sequence deduplication across samples, which dramatically reduces alignment time.

Workflow:
1. Preprocess all samples in parallel (trim, merge, collapse)
2. Collect unique sequences across all samples
3. Run triple-alignment ONCE on unique sequences
4. Apply classifications back to per-sample counts

Author: Kevin R. Roy
"""

import gzip
import logging
import tempfile
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

from ..integrations.crispresso import CRISPRessoResult, CRISPRessoRunner
from ..preprocessing.detection import run_auto_detection
from ..preprocessing.preprocess import run_preprocessing

logger = logging.getLogger(__name__)


def _classify_sequence_batch_worker(batch_data):
    """
    Worker function for parallel k-mer classification.

    This is a module-level function (not a method) for efficient pickling
    when using ProcessPoolExecutor.

    Args:
        batch_data: Tuple of (sequences, classifier, batch_idx)

    Returns:
        Dict mapping sequences to Classification objects
    """
    sequences, classifier, batch_idx = batch_data
    results = {}

    for seq in sequences:
        result = classifier.classify_sequence(seq)

        if result.is_wt:
            results[seq] = Classification(outcome='WT')
        elif result.best_template:
            results[seq] = Classification(
                outcome='HDR',
                template_id=result.best_template,
                is_perfect=True,  # Simplified for now
            )
        elif result.is_contamination:
            results[seq] = Classification(outcome='CONTAMINATION')
        elif result.is_ambiguous:
            results[seq] = Classification(outcome='AMBIGUOUS')
        else:
            results[seq] = Classification(outcome='UNKNOWN')

    return results


@dataclass
class SampleInfo:
    """Information about a sample for batch processing."""
    sample_id: str
    r1_path: Path
    r2_path: Optional[Path] = None
    expected_template: Optional[str] = None
    metadata: Dict[str, Any] = field(default_factory=dict)


@dataclass
class CollapsedSample:
    """Collapsed sequences from a single sample."""
    sample_id: str
    sequences: Dict[str, int]  # sequence → count
    total_reads: int

    @property
    def unique_sequences(self) -> int:
        """Number of unique sequences."""
        return len(self.sequences)


@dataclass
class Classification:
    """Classification result for a unique sequence."""
    outcome: str  # 'WT', 'HDR', 'NHEJ', 'LARGE_DELETION', 'AMBIGUOUS', 'UNKNOWN'
    template_id: Optional[str] = None  # For HDR, which template
    is_perfect: bool = False  # For HDR, perfect vs imperfect
    indel_info: Optional[str] = None  # For NHEJ, description of indel


@dataclass
class SampleResult:
    """Results for a single sample after batch processing."""
    sample_id: str
    total_reads: int
    unique_sequences: int

    # Outcome counts (k-mer + alignment)
    wt_count: int = 0
    hdr_count: int = 0
    hdr_perfect_count: int = 0
    hdr_imperfect_count: int = 0
    nhej_count: int = 0
    large_del_count: int = 0
    ambiguous_count: int = 0
    unknown_count: int = 0

    # Per-template HDR counts
    hdr_by_template: Dict[str, int] = field(default_factory=dict)

    # Expected template info
    expected_template: Optional[str] = None
    expected_template_count: int = 0

    # CRISPResso2 results (optional)
    crispresso_success: bool = False
    crispresso_total_reads: int = 0
    crispresso_aligned_reads: int = 0
    crispresso_wt_count: int = 0
    crispresso_hdr_count: int = 0
    crispresso_nhej_count: int = 0
    crispresso_error: Optional[str] = None

    @property
    def wt_rate(self) -> float:
        return self.wt_count / self.total_reads if self.total_reads > 0 else 0

    @property
    def hdr_rate(self) -> float:
        return self.hdr_count / self.total_reads if self.total_reads > 0 else 0

    @property
    def expected_template_rate(self) -> float:
        return self.expected_template_count / self.total_reads if self.total_reads > 0 else 0

    @property
    def crispresso_wt_rate(self) -> float:
        return self.crispresso_wt_count / self.crispresso_total_reads if self.crispresso_total_reads > 0 else 0

    @property
    def crispresso_hdr_rate(self) -> float:
        return self.crispresso_hdr_count / self.crispresso_total_reads if self.crispresso_total_reads > 0 else 0

    @property
    def crispresso_nhej_rate(self) -> float:
        return self.crispresso_nhej_count / self.crispresso_total_reads if self.crispresso_total_reads > 0 else 0


class BatchMultiTemplateProcessor:
    """
    Process multiple samples with global sequence deduplication.

    This processor is optimized for barcoded CRISPR experiments where:
    - Many samples share the same HDR templates
    - Same sequences appear across multiple samples
    - PCR duplicates are common within and across samples

    IMPORTANT: All samples must share the same reference sequence.
    For multi-locus experiments, create separate processor instances per reference.

    Workflow:
    1. Preprocess all samples in parallel (optional: trim, merge, collapse)
    2. Collect unique sequences across all samples
    3. Classify unique sequences using k-mer classification
    4. Apply classifications back to per-sample counts

    Example usage:
        processor = BatchMultiTemplateProcessor(
            reference='ATCG...',
            hdr_templates={'bc1': 'ATCG...bc1...', 'bc2': 'ATCG...bc2...'},
            guide='GAGTCCGAGCAGAAGAAGAA',
        )

        samples = [
            SampleInfo('sample_1', Path('s1_R1.fq.gz'), Path('s1_R2.fq.gz'), 'bc1'),
            SampleInfo('sample_2', Path('s2_R1.fq.gz'), Path('s2_R2.fq.gz'), 'bc2'),
        ]

        results = processor.process_batch(samples, output_dir=Path('results'))
    """

    def __init__(
        self,
        reference: str,
        hdr_templates: Dict[str, str],
        guide: str,
        kmer_size: int = 12,
        max_workers: int = 8,
        kmer_workers: int = 16,
        use_cache: bool = False,
        save_intermediates: bool = True,
    ):
        """
        Initialize processor for a single reference.

        Args:
            reference: Reference sequence (all samples must align to this)
            hdr_templates: Dict mapping template_id → HDR sequence
            guide: Guide RNA sequence (for NHEJ detection near cut site)
            kmer_size: Initial k-mer size (auto-increased for barcode templates)
            max_workers: Maximum parallel workers for preprocessing
            kmer_workers: Number of parallel workers for k-mer classification
            use_cache: Whether to use cached intermediate results if available
            save_intermediates: Whether to save intermediate outputs (k-mer results, collapsed sequences)
        """
        self.reference = reference
        self.hdr_templates = hdr_templates
        self.guide = guide
        self.kmer_size = kmer_size
        self.max_workers = max_workers
        self.kmer_workers = kmer_workers
        self.use_cache = use_cache
        self.save_intermediates = save_intermediates

        # Create k-mer classifier (lazy initialization)
        self._classifier = None

    @property
    def classifier(self):
        """Lazy initialization of k-mer classifier."""
        if self._classifier is None:
            from ..core.kmer import MultiTemplateKmerClassifier
            self._classifier = MultiTemplateKmerClassifier.from_multi_template(
                reference=self.reference,
                hdr_templates=self.hdr_templates,
                edit_positions_by_template={},  # Auto-detected for barcode templates
                kmer_size=self.kmer_size,
            )
            logger.info(f"Created k-mer classifier with k={self._classifier.kmer_size}")
        return self._classifier

    def process_batch(
        self,
        samples: List[SampleInfo],
        output_dir: Optional[Path] = None,
        enable_preprocessing: bool = False,
        amplicon_length: int = 200,
        max_reads: int = 100000,
        min_count: int = 1,
        max_unique_seqs: Optional[int] = None,
        min_global_count: int = 1,  # Filter by global count across all samples
        enable_alignment: bool = True,  # Alignment for NHEJ detection (default: True)
        alignment_only: bool = True,  # Pure alignment-based classification (default: True)
        enable_crispresso: bool = False,  # Optional CRISPResso2 comparison
    ) -> List[SampleResult]:
        """
        Process a batch of samples with global deduplication.

        Args:
            samples: List of SampleInfo objects
            output_dir: Optional output directory for results
            enable_preprocessing: If True, run auto-detection and preprocessing
                (UMI dedup, trim, merge, collapse) per sample
            amplicon_length: Expected amplicon length for auto-detection
            max_reads: Maximum reads to process per sample
            min_count: Minimum count threshold for unique sequences PER SAMPLE (default: 1).
                Sequences with fewer counts are filtered out to reduce noise.
            max_unique_seqs: Maximum unique sequences to keep per sample (default: None).
                If set, keeps only the top N sequences by count to speed up processing.
            min_global_count: Minimum total count across ALL samples (default: 1).
                Sequences appearing fewer than this many times globally are filtered out.
                Set to 2 to remove singletons (likely sequencing errors). This filter
                is applied AFTER building the global sequence map.
            enable_alignment: If True, run triple alignment on UNKNOWN sequences for
                NHEJ and large deletion detection (default: True). Only used when
                alignment_only=False (k-mer mode).
            alignment_only: Use pure alignment-based classification (default: True).
                Aligns all sequences to multi-reference FASTA (WT + all HDR variants)
                and classifies based on best alignment. Set to False for legacy k-mer
                mode (faster but less accurate for detecting imperfect HDR).

        Returns:
            List of SampleResult objects, one per sample
        """
        logger.info(f"Processing batch of {len(samples)} samples")
        if alignment_only:
            logger.info("Classification mode: ALIGNMENT-ONLY (multi-reference FASTA)")
        else:
            logger.info(f"Classification mode: HYBRID (k-mer + {'alignment' if enable_alignment else 'no alignment'})")
        if min_count > 1:
            logger.info(f"Filtering sequences with per-sample count < {min_count}")
        if max_unique_seqs:
            logger.info(f"Limiting to top {max_unique_seqs} unique sequences per sample")
        if min_global_count > 1:
            logger.info(f"Filtering sequences with global count < {min_global_count} (removes singletons/low-count errors)")

        # Phase 1: Preprocess and collapse sequences per sample
        if enable_preprocessing and output_dir:
            preprocess_dir = output_dir / "preprocessing"
            collapsed_samples = self._preprocess_and_collapse_samples(
                samples, preprocess_dir, amplicon_length, max_reads,
                min_count=min_count, max_unique_seqs=max_unique_seqs,
            )
        else:
            collapsed_samples = self._read_and_collapse_samples(
                samples, max_reads,
                min_count=min_count, max_unique_seqs=max_unique_seqs,
            )

        # Phase 2: Build global unique sequences mapping
        global_seqs, seq_to_samples = self._build_global_sequences(collapsed_samples)
        logger.info(f"Global unique sequences: {len(global_seqs)}")

        # Filter by global count (remove singletons/low-count errors)
        if min_global_count > 1:
            global_seqs, seq_to_samples = self._filter_global_sequences_by_count(
                global_seqs, seq_to_samples, min_global_count
            )

        # Save global sequences (if intermediate saving enabled)
        if output_dir:
            self._save_global_sequences(global_seqs, output_dir)

        # Phase 3: Classify unique sequences
        temp_dir = output_dir / "temp" if output_dir else Path(tempfile.mkdtemp())

        if alignment_only:
            # Pure alignment-based classification (multi-reference FASTA)
            classifications = self._classify_sequences_alignment_only(
                global_seqs,
                temp_dir=temp_dir,
                output_dir=output_dir,
            )
        else:
            # Hybrid k-mer + alignment classification
            classifications = self._classify_sequences(
                global_seqs,
                enable_alignment=enable_alignment,
                temp_dir=temp_dir,
                output_dir=output_dir,
            )

        # Phase 4: Apply classifications to per-sample counts
        results = self._apply_classifications(
            collapsed_samples, seq_to_samples, classifications, samples
        )

        # Phase 5: Optional CRISPResso2 comparison (per-sample on original FASTQs)
        if enable_crispresso and output_dir:
            logger.info("Running CRISPResso2 analysis per-sample...")
            crispresso_dir = output_dir / "crispresso"
            crispresso_results = self._run_crispresso_batch(samples, crispresso_dir)
            # Merge CRISPResso2 results into SampleResult objects
            self._merge_crispresso_results(results, crispresso_results)

        # Save results if output_dir provided
        if output_dir:
            self._save_results(results, output_dir)

        return results

    def _read_and_collapse_samples(
        self,
        samples: List[SampleInfo],
        max_reads: int = 100000,
        min_count: int = 1,
        max_unique_seqs: Optional[int] = None,
    ) -> List[CollapsedSample]:
        """
        Read FASTQs and collapse to unique sequences per sample.

        For now, this reads R1 only. In production, this would include
        adapter trimming and read merging.

        Args:
            samples: List of SampleInfo objects
            max_reads: Maximum reads to process per sample
            min_count: Minimum count threshold for unique sequences
            max_unique_seqs: Maximum unique sequences to keep per sample
        """
        collapsed = []

        for sample in samples:
            seq_counts: Dict[str, int] = Counter()
            total_reads = 0

            if not sample.r1_path.exists():
                logger.warning(f"FASTQ not found: {sample.r1_path}")
                continue

            opener = gzip.open if str(sample.r1_path).endswith('.gz') else open
            with opener(sample.r1_path, 'rt') as f:
                while total_reads < max_reads:
                    header = f.readline()
                    if not header:
                        break
                    seq = f.readline().strip().upper()
                    f.readline()  # +
                    f.readline()  # quality

                    total_reads += 1
                    seq_counts[seq] += 1

            # Apply sequence filters
            filtered_counts = self._apply_sequence_filters(
                seq_counts, min_count, max_unique_seqs
            )

            collapsed.append(CollapsedSample(
                sample_id=sample.sample_id,
                sequences=filtered_counts,
                total_reads=total_reads,
            ))

            logger.debug(
                f"{sample.sample_id}: {total_reads} reads → "
                f"{len(seq_counts)} unique → {len(filtered_counts)} filtered"
            )

        return collapsed

    def _preprocess_and_collapse_samples(
        self,
        samples: List[SampleInfo],
        preprocess_dir: Path,
        amplicon_length: int,
        max_reads: int,
        min_count: int = 1,
        max_unique_seqs: Optional[int] = None,
    ) -> List[CollapsedSample]:
        """
        Preprocess samples with auto-detection and collapse to unique sequences.

        For each sample:
        1. Auto-detect: library type, UMI presence, read overlap
        2. Run appropriate preprocessing:
           - TruSeq + UMI: dedup → trim → merge → collapse
           - TruSeq + no UMI: trim → merge → collapse
        3. Read collapsed/merged FASTQ and extract unique sequences

        Args:
            samples: List of SampleInfo objects
            preprocess_dir: Directory for preprocessing output
            amplicon_length: Expected amplicon length
            max_reads: Maximum reads to process
            min_count: Minimum count threshold for unique sequences
            max_unique_seqs: Maximum unique sequences to keep per sample

        Returns:
            List of CollapsedSample objects
        """
        preprocess_dir.mkdir(parents=True, exist_ok=True)
        collapsed = []

        for sample in samples:
            sample_dir = preprocess_dir / sample.sample_id
            sample_dir.mkdir(parents=True, exist_ok=True)

            if not sample.r1_path.exists():
                logger.warning(f"FASTQ not found: {sample.r1_path}")
                continue

            try:
                # Auto-detect library characteristics
                detection = run_auto_detection(
                    r1_path=sample.r1_path,
                    r2_path=sample.r2_path,
                    amplicon_length=amplicon_length,
                    reference=self.reference,
                )

                logger.info(
                    f"{sample.sample_id}: {detection.library.library_type}, "
                    f"UMI={detection.umi.has_umi if detection.umi else False}, "
                    f"Merge={detection.merge.should_merge if detection.merge else False}"
                )

                # Run preprocessing
                preprocess_result = run_preprocessing(
                    r1_path=sample.r1_path,
                    r2_path=sample.r2_path,
                    output_dir=sample_dir,
                    detection=detection,
                    threads=4,
                )

                if not preprocess_result.success:
                    logger.warning(
                        f"{sample.sample_id}: Preprocessing failed: "
                        f"{preprocess_result.error_message}"
                    )
                    continue

                # Read from preprocessed output
                if preprocess_result.output_format == 'merged' and preprocess_result.output_merged:
                    input_path = preprocess_result.output_merged
                else:
                    input_path = preprocess_result.output_r1

                # Collapse sequences from preprocessed output
                seq_counts, total_reads = self._collapse_fastq(input_path, max_reads)

                # Apply sequence filters
                filtered_counts = self._apply_sequence_filters(
                    seq_counts, min_count, max_unique_seqs
                )

                collapsed.append(CollapsedSample(
                    sample_id=sample.sample_id,
                    sequences=filtered_counts,
                    total_reads=total_reads,
                ))

                logger.info(
                    f"{sample.sample_id}: {preprocess_result.reads_before} raw → "
                    f"{preprocess_result.reads_after} preprocessed → "
                    f"{len(seq_counts)} unique → {len(filtered_counts)} filtered"
                )

            except Exception as e:
                logger.error(f"{sample.sample_id}: Error during preprocessing: {e}")
                # Fallback to raw FASTQ reading
                seq_counts, total_reads = self._collapse_fastq(sample.r1_path, max_reads)
                filtered_counts = self._apply_sequence_filters(
                    seq_counts, min_count, max_unique_seqs
                )
                collapsed.append(CollapsedSample(
                    sample_id=sample.sample_id,
                    sequences=filtered_counts,
                    total_reads=total_reads,
                ))

        return collapsed

    def _collapse_fastq(
        self,
        fastq_path: Path,
        max_reads: int,
    ) -> Tuple[Dict[str, int], int]:
        """
        Read FASTQ and collapse to unique sequences with counts.

        Supports pre-collapsed FASTQs with counts in headers (;count=N format).
        If no count is present in the header, counts occurrences of each sequence.

        Args:
            fastq_path: Path to FASTQ file
            max_reads: Maximum reads to process

        Returns:
            Tuple of (sequence → count dict, total reads)
        """
        import re
        count_pattern = re.compile(r';count=(\d+)')

        seq_counts: Dict[str, int] = Counter()
        total_reads = 0

        opener = gzip.open if str(fastq_path).endswith('.gz') else open
        with opener(fastq_path, 'rt') as f:
            while total_reads < max_reads:
                header = f.readline()
                if not header:
                    break
                seq = f.readline().strip().upper()
                f.readline()  # +
                f.readline()  # quality

                # Check if header contains a count (from pre-collapsed FASTQ)
                count_match = count_pattern.search(header)
                if count_match:
                    count = int(count_match.group(1))
                else:
                    count = 1

                total_reads += count
                seq_counts[seq] += count

        return dict(seq_counts), total_reads

    def _apply_sequence_filters(
        self,
        seq_counts: Dict[str, int],
        min_count: int = 1,
        max_unique_seqs: Optional[int] = None,
    ) -> Dict[str, int]:
        """
        Apply count and cardinality filters to sequence counts.

        Args:
            seq_counts: Dict mapping sequence → count
            min_count: Minimum count threshold (filter sequences below this)
            max_unique_seqs: Maximum unique sequences to keep (by count)

        Returns:
            Filtered dict mapping sequence → count
        """
        # Filter by minimum count
        if min_count > 1:
            seq_counts = {seq: count for seq, count in seq_counts.items() if count >= min_count}

        # Filter by max unique sequences (keep top N by count)
        if max_unique_seqs and len(seq_counts) > max_unique_seqs:
            # Sort by count descending and keep top N
            sorted_seqs = sorted(seq_counts.items(), key=lambda x: x[1], reverse=True)
            seq_counts = dict(sorted_seqs[:max_unique_seqs])

        return seq_counts

    def _build_global_sequences(
        self,
        collapsed_samples: List[CollapsedSample],
    ) -> Tuple[Set[str], Dict[str, Dict[str, int]]]:
        """
        Build global set of unique sequences and mapping to samples.

        Returns:
            global_seqs: Set of all unique sequences
            seq_to_samples: Dict mapping sequence → {sample_id: count}
        """
        global_seqs: Set[str] = set()
        seq_to_samples: Dict[str, Dict[str, int]] = {}

        for sample in collapsed_samples:
            for seq, count in sample.sequences.items():
                global_seqs.add(seq)
                if seq not in seq_to_samples:
                    seq_to_samples[seq] = {}
                seq_to_samples[seq][sample.sample_id] = count

        # Log statistics
        total_occurrences = sum(
            sum(counts.values())
            for counts in seq_to_samples.values()
        )
        logger.info(
            f"Built global sequence map: {len(global_seqs)} unique sequences, "
            f"{total_occurrences} total occurrences"
        )

        return global_seqs, seq_to_samples

    def _filter_global_sequences_by_count(
        self,
        global_seqs: Set[str],
        seq_to_samples: Dict[str, Dict[str, int]],
        min_global_count: int = 2,
    ) -> Tuple[Set[str], Dict[str, Dict[str, int]]]:
        """
        Filter out sequences with low global count across all samples.

        This is useful for removing likely sequencing errors (singletons) that
        appear only once across all samples.

        Args:
            global_seqs: Set of all unique sequences
            seq_to_samples: Dict mapping sequence → {sample_id: count}
            min_global_count: Minimum total count across all samples (default: 2)

        Returns:
            Filtered global_seqs and seq_to_samples
        """
        if min_global_count <= 1:
            return global_seqs, seq_to_samples

        # Calculate global count for each sequence
        seq_global_counts = {
            seq: sum(counts.values())
            for seq, counts in seq_to_samples.items()
        }

        # Filter sequences by global count
        filtered_seqs = {
            seq for seq, global_count in seq_global_counts.items()
            if global_count >= min_global_count
        }

        # Update seq_to_samples to only include filtered sequences
        filtered_seq_to_samples = {
            seq: counts for seq, counts in seq_to_samples.items()
            if seq in filtered_seqs
        }

        # Log statistics
        removed_count = len(global_seqs) - len(filtered_seqs)
        total_occurrences_before = sum(seq_global_counts.values())
        total_occurrences_after = sum(
            sum(counts.values())
            for counts in filtered_seq_to_samples.values()
        )
        removed_reads = total_occurrences_before - total_occurrences_after

        logger.info(
            f"Filtered global sequences by count >= {min_global_count}: "
            f"{len(global_seqs)} → {len(filtered_seqs)} unique sequences "
            f"({removed_count} removed, {removed_reads} reads discarded)"
        )

        return filtered_seqs, filtered_seq_to_samples

    def _save_kmer_results(
        self,
        classifications: Dict[str, Classification],
        output_dir: Path,
        filename: str = "kmer_classifications.pkl"
    ):
        """Save k-mer classification results to pickle for caching."""
        import pickle

        if not self.save_intermediates:
            return

        cache_file = output_dir / filename
        logger.info(f"Saving k-mer classifications to {cache_file}")

        with open(cache_file, 'wb') as f:
            pickle.dump(classifications, f)

        # Also save human-readable TSV
        tsv_file = output_dir / "kmer_classifications.tsv"
        with open(tsv_file, 'w') as f:
            f.write("sequence\toutcome\ttemplate_id\tis_perfect\tindel_info\n")
            for seq, result in classifications.items():
                f.write(f"{seq}\t{result.outcome}\t{result.template_id or ''}\t{result.is_perfect}\t{result.indel_info or ''}\n")

        logger.info(f"Saved {len(classifications)} classifications")

    def _load_kmer_results(
        self,
        output_dir: Path,
        filename: str = "kmer_classifications.pkl"
    ) -> Optional[Dict[str, Classification]]:
        """Load cached k-mer classification results if available."""
        import pickle

        if not self.use_cache:
            return None

        cache_file = output_dir / filename
        if not cache_file.exists():
            return None

        logger.info(f"Loading cached k-mer classifications from {cache_file}")
        with open(cache_file, 'rb') as f:
            classifications = pickle.load(f)

        logger.info(f"Loaded {len(classifications)} cached classifications")
        return classifications

    def _save_collapsed_sequences(
        self,
        sample_id: str,
        collapsed_seqs: Dict[str, int],
        output_dir: Path
    ):
        """Save collapsed sequences for a sample."""
        if not self.save_intermediates:
            return

        sample_dir = output_dir / "collapsed_sequences"
        sample_dir.mkdir(parents=True, exist_ok=True)

        output_file = sample_dir / f"{sample_id}_collapsed.tsv"

        with open(output_file, 'w') as f:
            f.write("sequence\tcount\n")
            for seq, count in sorted(collapsed_seqs.items(), key=lambda x: -x[1]):
                f.write(f"{seq}\t{count}\n")

    def _save_global_sequences(
        self,
        global_seqs: Set[str],
        output_dir: Path,
        filename: str = "global_unique_sequences.fasta"
    ):
        """Save global unique sequences to FASTA."""
        if not self.save_intermediates:
            return

        # Ensure output directory exists
        output_dir.mkdir(parents=True, exist_ok=True)

        output_file = output_dir / filename

        logger.info(f"Saving {len(global_seqs)} global unique sequences to {output_file}")

        with open(output_file, 'w') as f:
            for idx, seq in enumerate(sorted(global_seqs), 1):
                f.write(f">seq_{idx}\n{seq}\n")

    def _classify_sequences_parallel(
        self,
        sequences: Set[str],
        n_workers: Optional[int] = None
    ) -> Tuple[Dict[str, Classification], List[str]]:
        """
        Classify sequences in parallel using ProcessPoolExecutor.

        Args:
            sequences: Set of unique sequences to classify
            n_workers: Number of worker processes (default: self.kmer_workers)

        Returns:
            Tuple of (classifications dict, list of UNKNOWN sequences)
        """
        if n_workers is None:
            n_workers = self.kmer_workers

        logger.info(f"Running k-mer classification on {len(sequences)} unique sequences using {n_workers} workers...")

        # Convert set to list for batching
        seq_list = list(sequences)

        # Split sequences into batches
        batch_size = max(100, len(seq_list) // (n_workers * 4))  # ~4 batches per worker
        batches = [
            (seq_list[i:i+batch_size], self.classifier, batch_idx)
            for batch_idx, i in enumerate(range(0, len(seq_list), batch_size))
        ]

        logger.info(f"Split into {len(batches)} batches of ~{batch_size} sequences each")

        classifications = {}
        unknown_seqs = []

        # Process batches in parallel using module-level worker function
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            # Submit all batches (use module-level function for efficient pickling)
            future_to_batch = {
                executor.submit(_classify_sequence_batch_worker, batch): batch[2]
                for batch in batches
            }

            # Collect results as they complete
            for future in as_completed(future_to_batch):
                batch_idx = future_to_batch[future]
                try:
                    batch_results = future.result()
                    classifications.update(batch_results)

                    # Collect UNKNOWN sequences
                    for seq, classification in batch_results.items():
                        if classification.outcome == 'UNKNOWN':
                            unknown_seqs.append(seq)

                    logger.info(f"Completed batch {batch_idx + 1}/{len(batches)} ({len(classifications)}/{len(sequences)} total)")
                except Exception as e:
                    logger.error(f"Batch {batch_idx} failed: {e}")
                    raise

        # Log k-mer classification summary
        outcome_counts = Counter(c.outcome for c in classifications.values())
        logger.info(f"K-mer classification summary: {dict(outcome_counts)}")

        return classifications, unknown_seqs

    def _classify_sequences(
        self,
        sequences: Set[str],
        enable_alignment: bool = True,
        temp_dir: Optional[Path] = None,
        output_dir: Optional[Path] = None,
    ) -> Dict[str, Classification]:
        """
        Classify all unique sequences using hybrid k-mer + alignment approach.

        Workflow:
        1. Run k-mer classification on all sequences (fast, parallel)
        2. For UNKNOWN sequences, run triple alignment (complete NHEJ detection)
        3. Merge classifications

        This is the key optimization: classify each unique sequence ONCE,
        rather than classifying redundant copies across samples.

        Args:
            sequences: Set of unique sequences to classify
            enable_alignment: If True, run alignment on UNKNOWN sequences
            temp_dir: Temporary directory for alignment files
            output_dir: Output directory for saving intermediate results

        Returns:
            Dict mapping sequence → Classification
        """
        # Try to load cached k-mer results
        classifications = None
        unknown_seqs = []

        if output_dir and self.use_cache:
            classifications = self._load_kmer_results(output_dir)

        if classifications is None:
            # Phase 1: K-mer classification (parallel)
            classifications, unknown_seqs = self._classify_sequences_parallel(sequences)

            # Save k-mer results for future use
            if output_dir and self.save_intermediates:
                self._save_kmer_results(classifications, output_dir)
        else:
            # Extract UNKNOWN sequences from cached results
            logger.info("Using cached k-mer results, skipping classification")
            for seq, classification in classifications.items():
                if classification.outcome == 'UNKNOWN':
                    unknown_seqs.append(seq)

        # Phase 2: Alignment classification for UNKNOWN sequences (if enabled)
        if enable_alignment and unknown_seqs:
            logger.info(f"Running alignment classification on {len(unknown_seqs)} UNKNOWN sequences...")

            # Create alignment classifier
            from ..core.alignment_classifier import AlignmentClassifier, AlignmentConfig
            alignment_config = AlignmentConfig(threads=1)  # Safe for SLURM
            aligner = AlignmentClassifier(
                reference=self.reference,
                guide=self.guide,
                hdr_templates=self.hdr_templates,
                config=alignment_config,
            )

            # Classify UNKNOWN sequences
            alignment_dir = temp_dir / "alignment_classification" if temp_dir else Path(tempfile.mkdtemp())
            alignment_results = aligner.classify_sequences(set(unknown_seqs), alignment_dir)

            # Update classifications based on alignment results
            for seq, align_result in alignment_results.items():
                outcome = align_result.outcome.value  # Convert enum to string

                # Map alignment outcomes to Classification outcomes
                if 'nhej_insertion' in outcome:
                    classifications[seq] = Classification(outcome='NHEJ_INSERTION', indel_info=align_result.indel_info)
                elif 'nhej_deletion' in outcome:
                    classifications[seq] = Classification(outcome='NHEJ_DELETION', indel_info=align_result.indel_info)
                elif 'large_deletion' in outcome:
                    classifications[seq] = Classification(outcome='LARGE_DELETION', indel_info=align_result.indel_info)
                elif 'wild_type' in outcome:
                    classifications[seq] = Classification(outcome='WT')
                elif 'hdr' in outcome:
                    classifications[seq] = Classification(
                        outcome='HDR',
                        template_id=align_result.template_id,
                        is_perfect=align_result.is_perfect,
                    )
                elif 'unmapped' in outcome:
                    classifications[seq] = Classification(outcome='UNMAPPED')
                else:
                    # Keep as UNKNOWN
                    pass

            # Log alignment classification summary
            alignment_outcome_counts = Counter(c.outcome for c in classifications.values())
            logger.info(f"After alignment classification: {dict(alignment_outcome_counts)}")

        return classifications

    def _classify_sequences_alignment_only(
        self,
        sequences: Set[str],
        temp_dir: Path,
        output_dir: Optional[Path] = None,
    ) -> Dict[str, Classification]:
        """
        Classify sequences using PURE alignment-based approach.

        This replaces the hybrid k-mer + alignment approach with:
        1. Build multi-reference FASTA (WT + all HDR variants)
        2. Align all sequences to multi-ref using all 3 aligners
        3. Select best alignment per read
        4. Classify based on aligned reference + indel analysis

        Args:
            sequences: Set of unique sequences to classify
            temp_dir: Temporary directory for alignment files
            output_dir: Output directory for saving results

        Returns:
            Dict mapping sequence → Classification
        """
        from ..core.multi_ref_classifier import MultiRefClassifier
        from ..integrations.aligners import run_multi_ref_alignment
        from ..utils.multi_ref_builder import build_multi_reference_fasta

        logger.info(f"Running alignment-only classification on {len(sequences)} sequences")

        temp_dir.mkdir(parents=True, exist_ok=True)
        alignment_dir = temp_dir / "multi_ref_alignment"
        alignment_dir.mkdir(parents=True, exist_ok=True)

        # Step 1: Build multi-reference FASTA
        multi_ref_fasta = alignment_dir / "multi_reference.fasta"
        cut_sites = build_multi_reference_fasta(
            wt_reference=self.reference,
            hdr_templates=self.hdr_templates,
            output_path=multi_ref_fasta,
            guide_seq=self.guide,
        )

        logger.info(f"Built multi-reference FASTA with {len(self.hdr_templates) + 1} sequences")

        # Step 2: Write sequences to FASTQ
        query_fastq = alignment_dir / "query_sequences.fastq"
        with open(query_fastq, 'w') as f:
            for idx, seq in enumerate(sequences, 1):
                f.write(f"@seq_{idx}\n{seq}\n+\n{'I'*len(seq)}\n")

        # Step 3: Run alignment with all 3 aligners
        logger.info("Running triple alignment against multi-reference...")
        bam_paths = run_multi_ref_alignment(
            query_fastq=query_fastq,
            reference_fasta=multi_ref_fasta,
            output_dir=alignment_dir,
            threads=16,
        )

        logger.info(f"Alignment complete. BAMs: {list(bam_paths.keys())}")

        # Step 4: Classify using MultiRefClassifier
        classifier = MultiRefClassifier(
            cut_sites=cut_sites,
            cut_site_window=10,
            max_mismatches_imperfect=5,
            large_del_threshold=50,
        )

        align_classifications, summary = classifier.classify_from_bams(
            bam_paths=bam_paths,
            expected_barcode=None,  # Per-sample expected barcode handled later
        )

        logger.info(
            f"Classification summary: "
            f"WT={summary.wt_count}, HDR_perfect={summary.hdr_perfect_count}, "
            f"HDR_imperfect={summary.hdr_imperfect_count}, NHEJ={summary.nhej_count}, "
            f"large_del={summary.large_del_count}, unmapped={summary.unmapped_count}"
        )

        # Step 5: Convert to Classification objects
        classifications = {}

        # Build sequence index to map back from read names
        seq_list = list(sequences)
        seq_index = {f"seq_{idx+1}": seq for idx, seq in enumerate(seq_list)}

        for read_name, align_class in align_classifications.items():
            seq = seq_index.get(read_name)
            if seq is None:
                continue

            outcome = align_class.outcome
            if outcome == 'WT':
                classifications[seq] = Classification(outcome='WT')
            elif outcome == 'HDR_PERFECT':
                classifications[seq] = Classification(
                    outcome='HDR',
                    template_id=align_class.template_id,
                    is_perfect=True,
                )
            elif outcome == 'HDR_IMPERFECT':
                classifications[seq] = Classification(
                    outcome='HDR',
                    template_id=align_class.template_id,
                    is_perfect=False,
                )
            elif outcome == 'NHEJ':
                # Build indel info string
                indel_info = None
                if align_class.indels_near_cut:
                    indels = align_class.indels_near_cut
                    indel_info = "; ".join(
                        f"{'ins' if i.is_insertion else 'del'} {i.size}bp @ {i.position}"
                        for i in indels
                    )
                classifications[seq] = Classification(
                    outcome='NHEJ',
                    indel_info=indel_info,
                )
            elif outcome == 'LARGE_DEL':
                classifications[seq] = Classification(outcome='LARGE_DELETION')
            elif outcome == 'UNMAPPED':
                classifications[seq] = Classification(outcome='UNKNOWN')
            else:
                classifications[seq] = Classification(outcome='UNKNOWN')

        # Log final classification summary
        outcome_counts = Counter(c.outcome for c in classifications.values())
        logger.info(f"Final classification summary: {dict(outcome_counts)}")

        # Save classifications if output_dir provided
        if output_dir and self.save_intermediates:
            self._save_alignment_classifications(classifications, output_dir)

        return classifications

    def _save_alignment_classifications(
        self,
        classifications: Dict[str, Classification],
        output_dir: Path,
        filename: str = "alignment_classifications.tsv"
    ):
        """Save alignment-based classification results."""
        output_dir.mkdir(parents=True, exist_ok=True)

        output_file = output_dir / filename
        with open(output_file, 'w') as f:
            f.write("sequence\toutcome\ttemplate_id\tis_perfect\tindel_info\n")
            for seq, result in classifications.items():
                f.write(f"{seq}\t{result.outcome}\t{result.template_id or ''}\t{result.is_perfect}\t{result.indel_info or ''}\n")

        logger.info(f"Saved {len(classifications)} alignment classifications to {output_file}")

    def _apply_classifications(
        self,
        collapsed_samples: List[CollapsedSample],
        seq_to_samples: Dict[str, Dict[str, int]],
        classifications: Dict[str, Classification],
        sample_infos: List[SampleInfo],
    ) -> List[SampleResult]:
        """
        Apply pre-computed classifications to per-sample counts.
        """
        # Build sample_id → expected_template mapping
        expected_templates = {
            s.sample_id: s.expected_template
            for s in sample_infos
        }

        results = []
        for sample in collapsed_samples:
            result = SampleResult(
                sample_id=sample.sample_id,
                total_reads=sample.total_reads,
                unique_sequences=sample.unique_sequences,
                expected_template=expected_templates.get(sample.sample_id),
            )

            # Apply classifications
            for seq, count in sample.sequences.items():
                classification = classifications.get(seq)
                if classification is None:
                    result.unknown_count += count
                    continue

                outcome = classification.outcome
                if outcome == 'WT':
                    result.wt_count += count
                elif outcome == 'HDR':
                    result.hdr_count += count
                    tid = classification.template_id
                    if tid:
                        result.hdr_by_template[tid] = result.hdr_by_template.get(tid, 0) + count
                        if tid == result.expected_template:
                            result.expected_template_count += count
                    if classification.is_perfect:
                        result.hdr_perfect_count += count
                    else:
                        result.hdr_imperfect_count += count
                elif outcome in ('NHEJ', 'NHEJ_INSERTION', 'NHEJ_DELETION'):
                    # NHEJ includes both insertions and deletions
                    result.nhej_count += count
                elif outcome == 'LARGE_DELETION':
                    result.large_del_count += count
                elif outcome == 'AMBIGUOUS':
                    result.ambiguous_count += count
                elif outcome == 'UNMAPPED':
                    result.unknown_count += count
                else:
                    result.unknown_count += count

            results.append(result)

        return results

    def _run_crispresso_batch(
        self,
        samples: List[SampleInfo],
        output_dir: Path,
    ) -> Dict[str, CRISPRessoResult]:
        """
        Run CRISPResso2 per-sample on original R1/R2 FASTQs.

        Args:
            samples: List of SampleInfo objects
            output_dir: Base output directory for CRISPResso2 results

        Returns:
            Dict mapping sample_id to CRISPRessoResult
        """
        output_dir.mkdir(parents=True, exist_ok=True)
        runner = CRISPRessoRunner()

        if not runner.is_available():
            logger.warning("CRISPResso2 is not available, skipping")
            return {}

        results = {}

        for sample in samples:
            sample_output = output_dir / sample.sample_id
            logger.info(f"Running CRISPResso2 for {sample.sample_id}")

            # Build amplicon sequence - for now use reference
            # In production, may need to handle multiple templates
            amplicon_seq = self.reference

            result = runner.run(
                r1_fastq=sample.r1_path,
                r2_fastq=sample.r2_path if sample.r2_path else None,
                amplicon_seq=amplicon_seq,
                guide_seq=self.guide,
                output_dir=sample_output,
                name=sample.sample_id,
                quantification_window=10,
                min_identity=70,
            )

            results[sample.sample_id] = result

            if result.success:
                logger.info(
                    f"{sample.sample_id}: CRISPResso2 aligned {result.aligned_reads}/{result.total_reads} reads"
                )
            else:
                logger.warning(f"{sample.sample_id}: CRISPResso2 failed - {result.error_message}")

        return results

    def _merge_crispresso_results(
        self,
        results: List[SampleResult],
        crispresso_results: Dict[str, CRISPRessoResult],
    ):
        """
        Merge CRISPResso2 results into SampleResult objects (modifies in-place).

        Args:
            results: List of SampleResult objects
            crispresso_results: Dict mapping sample_id to CRISPRessoResult
        """
        for result in results:
            crispresso = crispresso_results.get(result.sample_id)
            if not crispresso:
                continue

            result.crispresso_success = crispresso.success
            result.crispresso_total_reads = crispresso.total_reads
            result.crispresso_aligned_reads = crispresso.aligned_reads

            if crispresso.success and crispresso.wt_rate is not None:
                # Convert rates back to counts
                result.crispresso_wt_count = int(crispresso.wt_rate * crispresso.total_reads)
                result.crispresso_hdr_count = int(crispresso.hdr_rate * crispresso.total_reads) if crispresso.hdr_rate else 0
                result.crispresso_nhej_count = int(crispresso.nhej_rate * crispresso.total_reads) if crispresso.nhej_rate else 0
            else:
                result.crispresso_error = crispresso.error_message

    def _save_results(
        self,
        results: List[SampleResult],
        output_dir: Path,
    ) -> None:
        """Save results to output directory."""
        output_dir.mkdir(parents=True, exist_ok=True)

        # Save summary table
        summary_path = output_dir / 'batch_summary.tsv'
        with open(summary_path, 'w') as f:
            headers = [
                'sample_id', 'total_reads', 'unique_sequences',
                'wt_count', 'wt_rate', 'hdr_count', 'hdr_rate',
                'nhej_count', 'unknown_count',
                'expected_template', 'expected_template_count', 'expected_template_rate',
                # CRISPResso2 results
                'crispresso_success', 'crispresso_total_reads', 'crispresso_aligned_reads',
                'crispresso_wt_count', 'crispresso_wt_rate',
                'crispresso_hdr_count', 'crispresso_hdr_rate',
                'crispresso_nhej_count', 'crispresso_nhej_rate',
                'crispresso_error',
            ]
            f.write('\t'.join(headers) + '\n')

            for r in results:
                row = [
                    r.sample_id, r.total_reads, r.unique_sequences,
                    r.wt_count, f'{r.wt_rate:.4f}',
                    r.hdr_count, f'{r.hdr_rate:.4f}',
                    r.nhej_count, r.unknown_count,
                    r.expected_template or '', r.expected_template_count,
                    f'{r.expected_template_rate:.4f}',
                    # CRISPResso2 results
                    r.crispresso_success, r.crispresso_total_reads, r.crispresso_aligned_reads,
                    r.crispresso_wt_count, f'{r.crispresso_wt_rate:.4f}',
                    r.crispresso_hdr_count, f'{r.crispresso_hdr_rate:.4f}',
                    r.crispresso_nhej_count, f'{r.crispresso_nhej_rate:.4f}',
                    r.crispresso_error or '',
                ]
                f.write('\t'.join(str(x) for x in row) + '\n')

        logger.info(f"Saved batch summary to {summary_path}")

        # Save per-template breakdown
        template_path = output_dir / 'per_template_hdr.tsv'
        with open(template_path, 'w') as f:
            f.write('sample_id\ttemplate_id\thdr_count\thdr_rate\tis_expected\n')
            for r in results:
                for tid, count in sorted(r.hdr_by_template.items()):
                    rate = count / r.total_reads if r.total_reads > 0 else 0
                    is_expected = tid == r.expected_template
                    f.write(f'{r.sample_id}\t{tid}\t{count}\t{rate:.4f}\t{is_expected}\n')

        logger.info(f"Saved per-template HDR to {template_path}")


def partition_samples_by_reference(
    samples: List[SampleInfo],
    sample_to_reference: Dict[str, str],
) -> Dict[str, List[SampleInfo]]:
    """
    Partition samples by their reference sequence.

    For multi-locus experiments, samples must be processed separately
    by reference to ensure valid alignments.

    Args:
        samples: List of all samples
        sample_to_reference: Dict mapping sample_id → reference_id

    Returns:
        Dict mapping reference_id → list of samples using that reference
    """
    partitions: Dict[str, List[SampleInfo]] = {}

    for sample in samples:
        ref_id = sample_to_reference.get(sample.sample_id, 'default')
        if ref_id not in partitions:
            partitions[ref_id] = []
        partitions[ref_id].append(sample)

    return partitions
