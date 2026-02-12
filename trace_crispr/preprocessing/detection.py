"""
Auto-detection of library type, read merging need, and CRISPResso mode.

Author: Kevin R. Roy
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Tuple, Optional
import gzip


# Adapter signatures
TN5_ADAPTER = "CTGTCTCTTATACACATCT"
TRUSEQ_ADAPTER = "AGATCGGAAGAGC"


@dataclass
class LibraryDetectionResult:
    """Result of library type detection."""
    library_type: str  # 'Tn5' or 'TruSeq'
    tn5_fraction: float
    truseq_fraction: float
    reads_checked: int
    confidence: str  # 'high', 'medium', 'low'

    def __str__(self) -> str:
        if self.library_type == 'Tn5':
            return f"Tn5 (Nextera adapter detected in {self.tn5_fraction:.1%} of reads)"
        else:
            return f"TruSeq (TruSeq adapter detected in {self.truseq_fraction:.1%} of reads)"


@dataclass
class MergeDetectionResult:
    """Result of merge need detection."""
    should_merge: bool
    estimated_overlap: float
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
    n_reads: int = 10000
) -> LibraryDetectionResult:
    """
    Auto-detect TruSeq vs Tn5 based on adapter signatures.

    Tn5 signature: CTGTCTCTTATACACATCT (Nextera/Tn5 transposase)
    TruSeq signature: AGATCGGAAGAGC

    Args:
        r1_path: Path to R1 FASTQ file
        r2_path: Optional path to R2 FASTQ file
        n_reads: Number of reads to sample

    Returns:
        LibraryDetectionResult with library type and confidence
    """
    tn5_count = 0
    truseq_count = 0
    reads_checked = 0

    open_func = gzip.open if str(r1_path).endswith('.gz') else open

    with open_func(r1_path, 'rt') as f:
        while reads_checked < n_reads:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()
            f.readline()  # +
            f.readline()  # qual

            if TN5_ADAPTER in seq:
                tn5_count += 1
            if TRUSEQ_ADAPTER in seq:
                truseq_count += 1

            reads_checked += 1

    if reads_checked == 0:
        return LibraryDetectionResult(
            library_type='TruSeq',
            tn5_fraction=0,
            truseq_fraction=0,
            reads_checked=0,
            confidence='low',
        )

    tn5_fraction = tn5_count / reads_checked
    truseq_fraction = truseq_count / reads_checked

    # Determine library type
    if tn5_fraction > truseq_fraction * 2:
        library_type = 'Tn5'
        confidence = 'high' if tn5_fraction > 0.3 else 'medium'
    elif truseq_fraction > tn5_fraction * 2:
        library_type = 'TruSeq'
        confidence = 'high' if truseq_fraction > 0.3 else 'medium'
    else:
        # Default to TruSeq if unclear
        library_type = 'TruSeq'
        confidence = 'low'

    return LibraryDetectionResult(
        library_type=library_type,
        tn5_fraction=tn5_fraction,
        truseq_fraction=truseq_fraction,
        reads_checked=reads_checked,
        confidence=confidence,
    )


def detect_merge_need(
    r1_path: Path,
    r2_path: Path,
    amplicon_length: int,
    n_reads: int = 1000
) -> MergeDetectionResult:
    """
    Detect if R1 and R2 should be merged based on overlap.

    Args:
        r1_path: Path to R1 FASTQ file
        r2_path: Path to R2 FASTQ file
        amplicon_length: Expected amplicon length
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
            estimated_overlap=0,
            reason="no reads found"
        )

    avg_read_length = sum(read_lengths) / len(read_lengths)

    # Simple heuristic: if reads could overlap by >=20bp, merge
    potential_overlap = (2 * avg_read_length) - amplicon_length

    if potential_overlap >= 20:
        overlap_fraction = potential_overlap / amplicon_length
        return MergeDetectionResult(
            should_merge=True,
            estimated_overlap=overlap_fraction,
            reason=f"estimated {overlap_fraction:.0%} overlap based on amplicon size"
        )
    else:
        return MergeDetectionResult(
            should_merge=False,
            estimated_overlap=max(0, potential_overlap / amplicon_length),
            reason="insufficient read overlap"
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
class AutoDetectionResults:
    """Combined auto-detection results."""
    library: LibraryDetectionResult
    merge: Optional[MergeDetectionResult]
    crispresso: CRISPRessoModeResult

    def print_summary(self):
        """Print auto-detection summary."""
        print("\nAuto-detection results:")
        print(f"  - Library type: {self.library}")
        if self.merge:
            print(f"  - Read merging: {self.merge}")
        print(f"  - CRISPResso mode: {self.crispresso}")


def run_auto_detection(
    r1_path: Path,
    r2_path: Optional[Path],
    amplicon_length: int,
) -> AutoDetectionResults:
    """
    Run all auto-detection steps.

    Args:
        r1_path: Path to R1 FASTQ
        r2_path: Optional path to R2 FASTQ
        amplicon_length: Expected amplicon length

    Returns:
        AutoDetectionResults with all detection results
    """
    # Detect library type
    library = detect_library_type(r1_path, r2_path)

    # Detect merge need if paired-end
    merge = None
    should_merge = False
    if r2_path:
        merge = detect_merge_need(r1_path, r2_path, amplicon_length)
        should_merge = merge.should_merge

    # Detect CRISPResso mode
    crispresso = detect_crispresso_mode(r1_path, r2_path, should_merge)

    return AutoDetectionResults(
        library=library,
        merge=merge,
        crispresso=crispresso,
    )
