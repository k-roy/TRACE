"""
Auto-detection of library type, read merging need, CRISPResso mode, and UMI length.

Author: Kevin R. Roy
"""

import gzip
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# Adapter signatures
TN5_ADAPTER = "CTGTCTCTTATACACATCT"
TRUSEQ_ADAPTER = "AGATCGGAAGAGC"


@dataclass
class LibraryDetectionResult:
    """Result of library type detection."""
    library_type: str  # 'Tn5' or 'TruSeq'
    reads_checked: int
    confidence: str  # 'high', 'medium', 'low'
    # Primary evidence: alignment position distribution
    position_clustering: float  # Fraction of reads at dominant start position (±5bp)
    n_unique_positions: int  # Number of distinct start positions observed
    # Secondary evidence: adapter signatures
    tn5_adapter_fraction: float
    truseq_adapter_fraction: float
    detection_reason: str  # Human-readable explanation

    def __str__(self) -> str:
        if self.library_type == 'Tn5':
            return f"Tn5 ({self.detection_reason})"
        else:
            return f"TruSeq ({self.detection_reason})"


@dataclass
class MergeDetectionResult:
    """Result of merge need detection."""
    should_merge: bool
    estimated_overlap_bp: int  # Estimated overlap in base pairs
    estimated_overlap_fraction: float  # Overlap as fraction of amplicon
    reason: str

    def __str__(self) -> str:
        if self.should_merge:
            return f"Enabled ({self.reason})"
        else:
            return f"Disabled ({self.reason})"


@dataclass
class CRISPRessoModeResult:
    """Result of CRISPResso mode detection."""
    mode: str  # 'single', 'paired', 'merged'
    reason: str

    def __str__(self) -> str:
        return f"{self.mode} ({self.reason})"


def detect_library_type(
    r1_path: Path,
    r2_path: Optional[Path] = None,
    reference: Optional[str] = None,
    n_reads: int = 1000
) -> LibraryDetectionResult:
    """
    Auto-detect TruSeq vs Tn5 based on alignment position distribution.

    Primary detection method:
    - TruSeq (fixed-target amplicon): Reads cluster at fixed start positions
      because primers bind at specific locations.
    - Tn5 (tagmented locus): Reads have scattered start positions because
      Tn5 transposase cuts randomly within the target region.

    Secondary evidence (adapter signatures):
    - Tn5: CTGTCTCTTATACACATCT (Nextera/Tn5 transposase)
    - TruSeq: AGATCGGAAGAGC

    Args:
        r1_path: Path to R1 FASTQ file
        r2_path: Optional path to R2 FASTQ file
        reference: Reference sequence for alignment-based detection
        n_reads: Number of reads to sample

    Returns:
        LibraryDetectionResult with library type and confidence
    """
    # Sample reads
    reads = []
    tn5_count = 0
    truseq_count = 0

    open_func = gzip.open if str(r1_path).endswith('.gz') else open

    with open_func(r1_path, 'rt') as f:
        while len(reads) < n_reads:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()
            f.readline()  # +
            f.readline()  # qual

            reads.append(seq)

            # Check for adapter signatures (secondary evidence)
            if TN5_ADAPTER in seq:
                tn5_count += 1
            if TRUSEQ_ADAPTER in seq:
                truseq_count += 1

    reads_checked = len(reads)
    if reads_checked == 0:
        return LibraryDetectionResult(
            library_type='TruSeq',
            reads_checked=0,
            confidence='low',
            position_clustering=0,
            n_unique_positions=0,
            tn5_adapter_fraction=0,
            truseq_adapter_fraction=0,
            detection_reason='no reads found, defaulting to TruSeq',
        )

    tn5_fraction = tn5_count / reads_checked
    truseq_fraction = truseq_count / reads_checked

    # Primary detection: alignment position distribution
    position_clustering = 0.0
    n_unique_positions = 0

    if reference:
        # Find alignment start positions using k-mer matching
        start_positions = _find_read_start_positions(reads, reference, kmer_size=15)

        if start_positions:
            position_counts = Counter(start_positions)
            n_unique_positions = len(position_counts)

            # Find the dominant position and count reads within ±5bp window
            if position_counts:
                dominant_pos, dominant_count = position_counts.most_common(1)[0]
                # Count reads within ±5bp of dominant position
                clustered_count = sum(
                    count for pos, count in position_counts.items()
                    if abs(pos - dominant_pos) <= 5
                )
                position_clustering = clustered_count / len(start_positions)

    # Decision logic:
    # Primary: Use position clustering if reference was provided
    # Secondary: Fall back to adapter signatures
    if reference and len(reads) > 100:
        # Position-based detection (primary)
        if position_clustering >= 0.7:
            # Most reads start at same position -> fixed primers -> TruSeq
            library_type = 'TruSeq'
            confidence = 'high' if position_clustering >= 0.85 else 'medium'
            detection_reason = f"{position_clustering:.0%} of reads cluster at fixed start position"
        elif n_unique_positions >= 20 and position_clustering < 0.5:
            # Many different start positions -> random tagmentation -> Tn5
            library_type = 'Tn5'
            confidence = 'high' if position_clustering < 0.3 else 'medium'
            detection_reason = f"scattered start positions ({n_unique_positions} unique positions)"
        else:
            # Ambiguous - use adapter signatures as tiebreaker
            if tn5_fraction > truseq_fraction * 2:
                library_type = 'Tn5'
                confidence = 'medium'
                detection_reason = f"Nextera adapter in {tn5_fraction:.1%} of reads"
            elif truseq_fraction > tn5_fraction * 2:
                library_type = 'TruSeq'
                confidence = 'medium'
                detection_reason = f"TruSeq adapter in {truseq_fraction:.1%} of reads"
            else:
                library_type = 'TruSeq'
                confidence = 'low'
                detection_reason = 'ambiguous, defaulting to TruSeq'
    else:
        # No reference provided - fall back to adapter-only detection
        if tn5_fraction > truseq_fraction * 2:
            library_type = 'Tn5'
            confidence = 'high' if tn5_fraction > 0.3 else 'medium'
            detection_reason = f"Nextera adapter in {tn5_fraction:.1%} of reads"
        elif truseq_fraction > tn5_fraction * 2:
            library_type = 'TruSeq'
            confidence = 'high' if truseq_fraction > 0.3 else 'medium'
            detection_reason = f"TruSeq adapter in {truseq_fraction:.1%} of reads"
        else:
            library_type = 'TruSeq'
            confidence = 'low'
            detection_reason = 'no adapters detected, defaulting to TruSeq'

    return LibraryDetectionResult(
        library_type=library_type,
        reads_checked=reads_checked,
        confidence=confidence,
        position_clustering=position_clustering,
        n_unique_positions=n_unique_positions,
        tn5_adapter_fraction=tn5_fraction,
        truseq_adapter_fraction=truseq_fraction,
        detection_reason=detection_reason,
    )


def _find_read_start_positions(
    reads: List[str],
    reference: str,
    kmer_size: int = 15
) -> List[int]:
    """
    Find where each read aligns to the reference using k-mer matching.

    Tries multiple k-mer positions in each read to handle UMI prefixes and
    other artifacts at read starts.

    Args:
        reads: List of read sequences
        reference: Reference sequence
        kmer_size: Size of k-mer to use for matching

    Returns:
        List of start positions (0-indexed) for reads that could be aligned
    """
    ref_upper = reference.upper()
    ref_rc = _reverse_complement(ref_upper)

    # Build k-mer index of reference (both strands)
    ref_kmers: Dict[str, List[Tuple[int, bool]]] = {}  # kmer -> [(pos, is_rc), ...]
    for i in range(len(ref_upper) - kmer_size + 1):
        kmer = ref_upper[i:i + kmer_size]
        if kmer not in ref_kmers:
            ref_kmers[kmer] = []
        ref_kmers[kmer].append((i, False))  # Forward strand

    # Also index reverse complement
    for i in range(len(ref_rc) - kmer_size + 1):
        kmer = ref_rc[i:i + kmer_size]
        if kmer not in ref_kmers:
            ref_kmers[kmer] = []
        # Position on original forward strand
        fwd_pos = len(ref_upper) - i - kmer_size
        ref_kmers[kmer].append((fwd_pos, True))  # RC strand

    # Find start positions for each read
    # Search comprehensively: check every position from 0-60 to handle both:
    # - TruSeq with UMIs (typically at fixed offsets like 0, 6, 8)
    # - Tn5 without UMIs (match can be at any offset due to adapter trimming)
    start_positions = []
    max_offset = 60

    for read in reads:
        if len(read) < kmer_size + 20:
            continue

        # Check every position from 0 to max_offset for comprehensive coverage
        for offset in range(0, min(max_offset, len(read) - kmer_size)):
            read_kmer = read[offset:offset + kmer_size].upper()
            if read_kmer in ref_kmers:
                # Found a match - calculate the read start position
                ref_pos, is_rc = ref_kmers[read_kmer][0]
                if is_rc:
                    # For RC alignment, the read start is at ref_pos + kmer_size - offset
                    # adjusted for reverse direction
                    read_start = ref_pos - offset
                else:
                    # For forward alignment, subtract the offset
                    read_start = ref_pos - offset

                if 0 <= read_start < len(ref_upper):
                    start_positions.append(read_start)
                    break

    return start_positions


def _reverse_complement(seq: str) -> str:
    """Return reverse complement of DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(b, 'N') for b in reversed(seq.upper()))


def detect_merge_need(
    r1_path: Path,
    r2_path: Path,
    amplicon_length: int,
    min_overlap_bp: int = 15,
    n_reads: int = 1000
) -> MergeDetectionResult:
    """
    Detect if R1 and R2 should be merged based on overlap.

    For TruSeq libraries, we use a threshold of 15bp overlap to enable merging.
    This allows better handling of reads that can be merged into a single
    higher-quality sequence.

    Args:
        r1_path: Path to R1 FASTQ file
        r2_path: Path to R2 FASTQ file
        amplicon_length: Expected amplicon length
        min_overlap_bp: Minimum overlap in bp to recommend merging (default: 15)
        n_reads: Number of reads to sample for length estimation

    Returns:
        MergeDetectionResult with merge recommendation
    """
    # Sample read lengths
    read_lengths = []
    open_func = gzip.open if str(r1_path).endswith('.gz') else open

    with open_func(r1_path, 'rt') as f:
        for i in range(min(n_reads, 1000)):
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()
            f.readline()  # +
            f.readline()  # qual
            read_lengths.append(len(seq))

    if not read_lengths:
        return MergeDetectionResult(
            should_merge=False,
            estimated_overlap_bp=0,
            estimated_overlap_fraction=0.0,
            reason="no reads found"
        )

    avg_read_length = sum(read_lengths) / len(read_lengths)

    # Calculate potential overlap: 2*read_length - amplicon_length
    potential_overlap_bp = int((2 * avg_read_length) - amplicon_length)
    overlap_fraction = potential_overlap_bp / amplicon_length if amplicon_length > 0 else 0

    if potential_overlap_bp >= min_overlap_bp:
        return MergeDetectionResult(
            should_merge=True,
            estimated_overlap_bp=potential_overlap_bp,
            estimated_overlap_fraction=overlap_fraction,
            reason=f"~{potential_overlap_bp}bp overlap ({overlap_fraction:.0%} of amplicon)"
        )
    else:
        return MergeDetectionResult(
            should_merge=False,
            estimated_overlap_bp=max(0, potential_overlap_bp),
            estimated_overlap_fraction=max(0, overlap_fraction),
            reason=f"insufficient overlap ({max(0, potential_overlap_bp)}bp < {min_overlap_bp}bp threshold)"
        )


def detect_crispresso_mode(
    r1_path: Path,
    r2_path: Optional[Path],
    should_merge: bool
) -> CRISPRessoModeResult:
    """
    Determine best CRISPResso2 mode based on data.

    Args:
        r1_path: Path to R1 FASTQ
        r2_path: Optional path to R2 FASTQ
        should_merge: Whether reads should be merged

    Returns:
        CRISPRessoModeResult with recommended mode
    """
    if r2_path is None:
        return CRISPRessoModeResult(
            mode='single',
            reason='single-end data'
        )

    if should_merge:
        return CRISPRessoModeResult(
            mode='merged',
            reason='overlapping paired-end reads'
        )
    else:
        return CRISPRessoModeResult(
            mode='paired',
            reason='non-overlapping paired-end reads'
        )


@dataclass
class UMIDetectionResult:
    """Result of UMI detection."""
    has_umi: bool
    r1_umi_length: int
    r2_umi_length: int
    r1_primer_seq: str  # Detected primer after UMI
    r2_primer_seq: str  # Detected primer after UMI
    confidence: str  # 'high', 'medium', 'low'
    reads_checked: int

    def __str__(self) -> str:
        if self.has_umi:
            if self.r1_umi_length == self.r2_umi_length:
                return f"UMIs of length {self.r1_umi_length} bp detected"
            else:
                return f"UMIs detected: R1={self.r1_umi_length} bp, R2={self.r2_umi_length} bp"
        else:
            return "No UMIs detected"


@dataclass
class PreprocessingMode:
    """Preprocessing protocol to apply based on detection results."""
    # Steps to perform (in order)
    dedup: bool  # UMI-based deduplication (TruSeq) or position-based (Tn5, post-align)
    trim: bool  # Adapter trimming
    merge: bool  # Read merging (requires overlap)
    collapse: bool  # Collapse identical sequences (post-merge)
    # Output format
    output_format: str  # 'merged' (single FASTQ) or 'paired' (R1+R2)
    description: str  # Human-readable description

    def __str__(self) -> str:
        steps = []
        if self.dedup:
            steps.append("dedup")
        if self.trim:
            steps.append("trim")
        if self.merge:
            steps.append("merge")
        if self.collapse:
            steps.append("collapse")
        return f"{'-'.join(steps)} → {self.output_format}"


@dataclass
class AutoDetectionResults:
    """Combined auto-detection results."""
    library: LibraryDetectionResult
    merge: Optional[MergeDetectionResult]
    crispresso: CRISPRessoModeResult
    umi: Optional[UMIDetectionResult] = None
    preprocessing: Optional[PreprocessingMode] = None

    def print_summary(self):
        """Print auto-detection summary."""
        print("\nAuto-detection results:")
        print(f"  - Library type: {self.library}")
        if self.umi and self.umi.has_umi:
            print(f"  - UMI detection: {self.umi}")
        if self.merge:
            print(f"  - Read overlap: {self.merge}")
        if self.preprocessing:
            print(f"  - Preprocessing: {self.preprocessing}")
            print(f"    ({self.preprocessing.description})")
        print(f"  - CRISPResso mode: {self.crispresso}")


def detect_umi_length(
    r1_path: Path,
    r2_path: Optional[Path],
    reference: str,
    n_reads: int = 1000
) -> UMIDetectionResult:
    """
    Auto-detect UMI length by finding the transition from high to low diversity.

    For TruSeq libraries with UMIs, the structure is:
    - [UMI (N bp)] [Gene-specific primer] [Amplicon sequence]

    UMIs are random sequences with high diversity, while the primer region
    has low diversity (same sequence in most reads). This function finds
    the position where diversity drops sharply.

    Args:
        r1_path: Path to R1 FASTQ file
        r2_path: Optional path to R2 FASTQ file
        reference: Reference amplicon sequence
        n_reads: Number of reads to sample

    Returns:
        UMIDetectionResult with detected UMI lengths
    """
    def find_umi_boundary_by_diversity(reads: List[str], max_umi_len: int = 15, kmer_size: int = 6) -> Tuple[int, str]:
        """Find UMI length by detecting where diversity drops.

        UMI region: high diversity (many unique k-mers)
        Primer region: low diversity (few unique k-mers, one dominates)

        Returns (umi_length, detected_primer_sequence)
        """
        if not reads:
            return 0, ""

        # Analyze diversity at each position
        diversity_scores = []
        top_kmers = []

        for pos in range(max_umi_len + 1):
            kmers = Counter()
            for read in reads:
                if pos + kmer_size <= len(read):
                    kmers[read[pos:pos + kmer_size].upper()] += 1

            if kmers:
                # Calculate diversity as fraction of unique k-mers
                n_unique = len(kmers)
                n_reads_with_kmer = sum(kmers.values())
                diversity = n_unique / n_reads_with_kmer if n_reads_with_kmer > 0 else 0

                # Get top k-mer
                top_kmer, top_count = kmers.most_common(1)[0]
                top_fraction = top_count / n_reads_with_kmer if n_reads_with_kmer > 0 else 0

                diversity_scores.append((pos, diversity, top_fraction))
                top_kmers.append((pos, top_kmer, top_count))
            else:
                diversity_scores.append((pos, 1.0, 0.0))
                top_kmers.append((pos, "", 0))

        # Find the position where diversity drops significantly (primer starts)
        # Use three-pass approach to handle variable primer quality
        # Pass 1: Strong signal (>50% consensus) - high confidence
        # Pass 2: Weak signal (>30% consensus, ≥4bp UMI) - medium confidence
        # Pass 3: Best available signal (≥5bp position) - low confidence fallback
        umi_length = 0

        # Pass 1: Strong threshold (0.5)
        for pos, diversity, top_fraction in diversity_scores:
            if top_fraction > 0.5 and pos > 0:
                umi_length = pos
                break

        # Pass 2: Weak threshold (0.3), require at least 4bp UMI
        if umi_length == 0:
            for pos, diversity, top_fraction in diversity_scores:
                if top_fraction > 0.3 and pos > 3:
                    umi_length = pos
                    break

        # Pass 3: Best available signal (for poor quality data)
        # Find FIRST position in range 4-8 where top_fraction shows significant increase
        if umi_length == 0:
            prev_frac = 0
            for pos, diversity, top_fraction in diversity_scores:
                if 4 <= pos <= 8:
                    # Look for jump of >10% indicating primer start
                    if top_fraction > 0.15 and (top_fraction - prev_frac) > 0.10:
                        umi_length = pos
                        break
                    prev_frac = top_fraction

            # If still no detection, default to 6bp (safe fallback for failed detection)
            if umi_length == 0 and len(diversity_scores) > 6:
                umi_length = 6

        # Get the primer sequence at the detected boundary
        primer_seq = ""
        if umi_length > 0 and umi_length < len(top_kmers):
            # Extract longer primer from reads at that position
            primer_kmers = Counter()
            primer_len = 15  # Look for 15-mer
            for read in reads:
                if umi_length + primer_len <= len(read):
                    primer_kmers[read[umi_length:umi_length + primer_len].upper()] += 1
            if primer_kmers:
                primer_seq = primer_kmers.most_common(1)[0][0]

        return umi_length, primer_seq

    def reverse_complement(seq: str) -> str:
        """Return reverse complement of DNA sequence."""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        return ''.join(complement.get(b, 'N') for b in reversed(seq.upper()))

    # Sample R1 reads
    r1_reads = []
    open_func = gzip.open if str(r1_path).endswith('.gz') else open
    with open_func(r1_path, 'rt') as f:
        count = 0
        while count < n_reads:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()
            f.readline()  # +
            f.readline()  # qual
            r1_reads.append(seq)
            count += 1

    r1_umi_length, r1_primer = find_umi_boundary_by_diversity(r1_reads)

    # Sample R2 reads if available
    r2_umi_length = 0
    r2_primer = ""
    if r2_path:
        r2_reads = []
        open_func = gzip.open if str(r2_path).endswith('.gz') else open
        with open_func(r2_path, 'rt') as f:
            count = 0
            while count < n_reads:
                header = f.readline().strip()
                if not header:
                    break
                seq = f.readline().strip()
                f.readline()  # +
                f.readline()  # qual
                r2_reads.append(seq)
                count += 1

        r2_umi_length, r2_primer = find_umi_boundary_by_diversity(r2_reads)

    # Determine if UMIs are present
    # UMI is typically 4-12 bp; consider it present if offset is >= 4
    has_umi = r1_umi_length >= 4 or r2_umi_length >= 4

    # Assess confidence
    if has_umi:
        if r1_umi_length == r2_umi_length and r1_umi_length >= 4:
            confidence = 'high'
        elif r1_umi_length >= 4 and r2_umi_length >= 4:
            confidence = 'medium'
        else:
            confidence = 'low'
    else:
        confidence = 'high'  # Confident no UMI

    return UMIDetectionResult(
        has_umi=has_umi,
        r1_umi_length=r1_umi_length,
        r2_umi_length=r2_umi_length,
        r1_primer_seq=r1_primer,
        r2_primer_seq=r2_primer,
        confidence=confidence,
        reads_checked=len(r1_reads),
    )


def _determine_preprocessing_mode(
    library_type: str,
    has_umi: bool,
    should_merge: bool,
    is_paired: bool,
) -> PreprocessingMode:
    """
    Determine preprocessing protocol based on detection results.

    TruSeq workflows:
    - UMI + overlap: dedup → trim → merge → collapse → merged FASTQ
    - UMI + no overlap: dedup → trim → paired FASTQ
    - no UMI + overlap: trim → merge → collapse → merged FASTQ
    - no UMI + no overlap: trim → paired FASTQ

    Tn5 workflows:
    - Always: trim → (paired or merged based on overlap) → align → position-based dedup
    """
    if library_type == 'TruSeq':
        if has_umi and should_merge:
            return PreprocessingMode(
                dedup=True,
                trim=True,
                merge=True,
                collapse=True,
                output_format='merged',
                description="UMI dedup → trim → merge → collapse"
            )
        elif has_umi and not should_merge:
            return PreprocessingMode(
                dedup=True,
                trim=True,
                merge=False,
                collapse=False,
                output_format='paired',
                description="UMI dedup → trim → paired-end mapping"
            )
        elif not has_umi and should_merge:
            return PreprocessingMode(
                dedup=False,
                trim=True,
                merge=True,
                collapse=True,
                output_format='merged',
                description="trim → merge → collapse"
            )
        else:  # no UMI, no overlap
            return PreprocessingMode(
                dedup=False,
                trim=True,
                merge=False,
                collapse=False,
                output_format='paired',
                description="trim → paired-end mapping"
            )
    else:  # Tn5
        if should_merge:
            return PreprocessingMode(
                dedup=False,  # Position-based dedup happens post-alignment for Tn5
                trim=True,
                merge=True,
                collapse=True,
                output_format='merged',
                description="trim → merge → collapse → align → position dedup"
            )
        elif is_paired:
            return PreprocessingMode(
                dedup=False,
                trim=True,
                merge=False,
                collapse=False,
                output_format='paired',
                description="trim → paired-end mapping → position dedup"
            )
        else:  # Single-end Tn5
            return PreprocessingMode(
                dedup=False,
                trim=True,
                merge=False,
                collapse=False,
                output_format='single',
                description="trim → single-end mapping → position dedup"
            )


def run_auto_detection(
    r1_path: Path,
    r2_path: Optional[Path],
    amplicon_length: int,
    reference: Optional[str] = None,
) -> AutoDetectionResults:
    """
    Run all auto-detection steps.

    Detection order for TruSeq:
    1. Library type (TruSeq vs Tn5)
    2. UMI presence and length (TruSeq only)
    3. Read overlap (>=15bp threshold)
    4. Determine preprocessing protocol

    Args:
        r1_path: Path to R1 FASTQ
        r2_path: Optional path to R2 FASTQ
        amplicon_length: Expected amplicon length
        reference: Optional reference sequence for position-based detection

    Returns:
        AutoDetectionResults with all detection results
    """
    # Step 1: Detect library type (use reference for position-based detection)
    library = detect_library_type(r1_path, r2_path, reference=reference)

    # Step 2: Detect UMI if TruSeq library and reference provided
    umi = None
    has_umi = False
    if library.library_type == 'TruSeq' and reference:
        umi = detect_umi_length(r1_path, r2_path, reference)
        has_umi = umi.has_umi

    # Step 3: Detect merge need if paired-end (using 15bp threshold)
    merge = None
    should_merge = False
    if r2_path:
        merge = detect_merge_need(r1_path, r2_path, amplicon_length, min_overlap_bp=15)
        should_merge = merge.should_merge

    # Step 4: Determine preprocessing protocol
    preprocessing = _determine_preprocessing_mode(
        library_type=library.library_type,
        has_umi=has_umi,
        should_merge=should_merge,
        is_paired=r2_path is not None,
    )

    # Detect CRISPResso mode based on preprocessing output
    if preprocessing.output_format == 'merged':
        crispresso = CRISPRessoModeResult(mode='merged', reason='merged reads from preprocessing')
    elif preprocessing.output_format == 'paired':
        crispresso = CRISPRessoModeResult(mode='paired', reason='paired-end from preprocessing')
    else:
        crispresso = CRISPRessoModeResult(mode='single', reason='single-end data')

    return AutoDetectionResults(
        library=library,
        merge=merge,
        crispresso=crispresso,
        umi=umi,
        preprocessing=preprocessing,
    )
