"""
Preprocessing orchestration for TRACE pipeline.

Handles different preprocessing modes based on library type, UMI presence, and overlap:
- TruSeq + UMI + overlap: dedup → trim → merge → collapse → merged FASTQ
- TruSeq + UMI + no overlap: dedup → trim → paired FASTQs
- TruSeq + no UMI + overlap: trim → merge → collapse → merged FASTQ
- TruSeq + no UMI + no overlap: trim → paired FASTQs
- Tn5: trim → (merged or paired) → align → position-based dedup (post-alignment)

Author: Kevin R. Roy
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple, Dict
import subprocess
import gzip
import logging
from collections import Counter

from .detection import PreprocessingMode, AutoDetectionResults

logger = logging.getLogger(__name__)


@dataclass
class PreprocessingResult:
    """Result of preprocessing pipeline."""
    success: bool
    output_r1: Optional[Path]  # For paired output
    output_r2: Optional[Path]  # For paired output (None if merged)
    output_merged: Optional[Path]  # For merged output
    output_format: str  # 'paired' or 'merged'
    reads_before: int
    reads_after: int
    duplicates_removed: int
    error_message: str = ""

    @property
    def dedup_rate(self) -> float:
        """Calculate deduplication rate."""
        if self.reads_before == 0:
            return 0.0
        return self.duplicates_removed / self.reads_before


def run_preprocessing(
    r1_path: Path,
    r2_path: Optional[Path],
    output_dir: Path,
    detection: AutoDetectionResults,
    threads: int = 4,
) -> PreprocessingResult:
    """
    Run preprocessing based on detection results.

    Args:
        r1_path: Path to R1 FASTQ
        r2_path: Path to R2 FASTQ (None for single-end)
        output_dir: Output directory for intermediate and final files
        detection: Auto-detection results with preprocessing mode
        threads: Number of threads for tools

    Returns:
        PreprocessingResult with output paths and statistics
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    mode = detection.preprocessing

    if mode is None:
        logger.warning("No preprocessing mode detected, using default (trim only)")
        return _run_trim_only(r1_path, r2_path, output_dir, threads)

    logger.info(f"Running preprocessing: {mode.description}")

    # Track statistics
    reads_before = count_reads(r1_path)
    duplicates_removed = 0

    # Current input files (will be updated as we go through steps)
    current_r1 = r1_path
    current_r2 = r2_path

    # Step 1: UMI-based deduplication (if applicable)
    if mode.dedup and detection.library.library_type == 'TruSeq' and detection.umi:
        logger.info("  Step 1: UMI-based deduplication...")
        dedup_dir = output_dir / "dedup"
        dedup_result = run_umi_dedup(
            current_r1, current_r2,
            dedup_dir,
            umi_length_r1=detection.umi.r1_umi_length,
            umi_length_r2=detection.umi.r2_umi_length if detection.umi else 0,
        )
        if dedup_result['success']:
            current_r1 = dedup_result['output_r1']
            current_r2 = dedup_result['output_r2']
            duplicates_removed = dedup_result['duplicates_removed']
        else:
            logger.warning(f"  Deduplication failed: {dedup_result.get('error', 'unknown')}")

    # Step 2: Adapter trimming
    if mode.trim:
        logger.info("  Step 2: Adapter trimming...")
        trim_dir = output_dir / "trimmed"
        from .trimming import trim_adapters
        trim_result = trim_adapters(
            current_r1, current_r2,
            trim_dir,
            library_type=detection.library.library_type.lower(),
            threads=threads,
        )
        if trim_result.success:
            current_r1 = trim_result.output_r1
            current_r2 = trim_result.output_r2
        else:
            logger.warning(f"  Trimming failed: {trim_result.error_message}")

    # Step 3: Read merging (if applicable)
    merged_path = None
    if mode.merge and current_r2:
        logger.info("  Step 3: Read merging...")
        merge_dir = output_dir / "merged"
        merge_result = run_read_merge(
            current_r1, current_r2,
            merge_dir,
            threads=threads,
        )
        if merge_result['success']:
            merged_path = merge_result['merged_path']
        else:
            logger.warning(f"  Merging failed: {merge_result.get('error', 'unknown')}")
            # Fall back to paired-end
            mode = PreprocessingMode(
                dedup=mode.dedup,
                trim=mode.trim,
                merge=False,
                collapse=False,
                output_format='paired',
                description=mode.description + " (merge failed, using paired)"
            )

    # Step 4: Collapse identical sequences (if applicable)
    if mode.collapse and merged_path:
        logger.info("  Step 4: Collapsing identical sequences...")
        collapse_dir = output_dir / "collapsed"
        collapse_result = run_collapse(merged_path, collapse_dir)
        if collapse_result['success']:
            merged_path = collapse_result['collapsed_path']
        else:
            logger.warning(f"  Collapse failed: {collapse_result.get('error', 'unknown')}")

    # Determine output
    reads_after = count_reads(merged_path if merged_path else current_r1)

    if mode.output_format == 'merged' and merged_path:
        return PreprocessingResult(
            success=True,
            output_r1=None,
            output_r2=None,
            output_merged=merged_path,
            output_format='merged',
            reads_before=reads_before,
            reads_after=reads_after,
            duplicates_removed=duplicates_removed,
        )
    else:
        return PreprocessingResult(
            success=True,
            output_r1=current_r1,
            output_r2=current_r2,
            output_merged=None,
            output_format='paired' if current_r2 else 'single',
            reads_before=reads_before,
            reads_after=reads_after,
            duplicates_removed=duplicates_removed,
        )


def _run_trim_only(
    r1_path: Path,
    r2_path: Optional[Path],
    output_dir: Path,
    threads: int,
) -> PreprocessingResult:
    """Run trimming only (fallback mode)."""
    from .trimming import trim_adapters

    reads_before = count_reads(r1_path)

    trim_result = trim_adapters(
        r1_path, r2_path,
        output_dir / "trimmed",
        library_type='auto',
        threads=threads,
    )

    if trim_result.success:
        reads_after = count_reads(trim_result.output_r1)
        return PreprocessingResult(
            success=True,
            output_r1=trim_result.output_r1,
            output_r2=trim_result.output_r2,
            output_merged=None,
            output_format='paired' if r2_path else 'single',
            reads_before=reads_before,
            reads_after=reads_after,
            duplicates_removed=0,
        )
    else:
        return PreprocessingResult(
            success=False,
            output_r1=r1_path,
            output_r2=r2_path,
            output_merged=None,
            output_format='paired' if r2_path else 'single',
            reads_before=reads_before,
            reads_after=reads_before,
            duplicates_removed=0,
            error_message=trim_result.error_message,
        )


def run_umi_dedup(
    r1_path: Path,
    r2_path: Optional[Path],
    output_dir: Path,
    umi_length_r1: int,
    umi_length_r2: int = 0,
) -> Dict:
    """
    Perform UMI-based deduplication.

    Extracts UMIs from the start of reads, uses them to identify duplicates,
    and keeps one representative read per UMI.

    Args:
        r1_path: Path to R1 FASTQ
        r2_path: Path to R2 FASTQ
        output_dir: Output directory
        umi_length_r1: UMI length at start of R1
        umi_length_r2: UMI length at start of R2

    Returns:
        Dict with 'success', 'output_r1', 'output_r2', 'duplicates_removed'
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    output_r1 = output_dir / f"{r1_path.stem.replace('.fastq', '')}_dedup.fastq.gz"
    output_r2 = output_dir / f"{r2_path.stem.replace('.fastq', '')}_dedup.fastq.gz" if r2_path else None

    # UMI -> (best_quality_read, qual_sum)
    umi_best: Dict[str, Tuple] = {}
    total_reads = 0

    open_r1 = gzip.open if str(r1_path).endswith('.gz') else open
    open_r2 = gzip.open if r2_path and str(r2_path).endswith('.gz') else open

    try:
        # First pass: find best read for each UMI
        with open_r1(r1_path, 'rt') as f1:
            f2 = open_r2(r2_path, 'rt') if r2_path else None

            while True:
                # Read R1
                header1 = f1.readline().strip()
                if not header1:
                    break
                seq1 = f1.readline().strip()
                plus1 = f1.readline().strip()
                qual1 = f1.readline().strip()

                # Read R2 if available
                if f2:
                    header2 = f2.readline().strip()
                    seq2 = f2.readline().strip()
                    plus2 = f2.readline().strip()
                    qual2 = f2.readline().strip()
                else:
                    header2 = seq2 = plus2 = qual2 = ""

                total_reads += 1

                # Extract UMIs
                umi_r1 = seq1[:umi_length_r1] if umi_length_r1 > 0 else ""
                umi_r2 = seq2[:umi_length_r2] if umi_length_r2 > 0 and r2_path else ""
                umi = umi_r1 + umi_r2

                # Trim UMIs from sequences
                trimmed_seq1 = seq1[umi_length_r1:]
                trimmed_qual1 = qual1[umi_length_r1:]
                trimmed_seq2 = seq2[umi_length_r2:] if r2_path else ""
                trimmed_qual2 = qual2[umi_length_r2:] if r2_path else ""

                # Calculate quality score
                qual_sum = sum(ord(c) - 33 for c in trimmed_qual1)
                if r2_path:
                    qual_sum += sum(ord(c) - 33 for c in trimmed_qual2)

                # Keep best quality read for each UMI
                if umi not in umi_best or qual_sum > umi_best[umi][1]:
                    umi_best[umi] = (
                        (header1, trimmed_seq1, trimmed_qual1,
                         header2, trimmed_seq2, trimmed_qual2),
                        qual_sum
                    )

            if f2:
                f2.close()

        # Write deduplicated reads
        with gzip.open(output_r1, 'wt') as out1:
            out2 = gzip.open(output_r2, 'wt') if output_r2 else None

            for umi, (read_data, _) in umi_best.items():
                header1, seq1, qual1, header2, seq2, qual2 = read_data

                out1.write(f"{header1}\n{seq1}\n+\n{qual1}\n")

                if out2 and header2:
                    out2.write(f"{header2}\n{seq2}\n+\n{qual2}\n")

            if out2:
                out2.close()

        duplicates = total_reads - len(umi_best)

        return {
            'success': True,
            'output_r1': output_r1,
            'output_r2': output_r2,
            'total_reads': total_reads,
            'unique_reads': len(umi_best),
            'duplicates_removed': duplicates,
        }

    except Exception as e:
        logger.error(f"UMI deduplication failed: {e}")
        return {
            'success': False,
            'error': str(e),
            'output_r1': r1_path,
            'output_r2': r2_path,
            'duplicates_removed': 0,
        }


def run_read_merge(
    r1_path: Path,
    r2_path: Path,
    output_dir: Path,
    threads: int = 4,
) -> Dict:
    """
    Merge overlapping paired-end reads using BBMerge.

    Args:
        r1_path: Path to R1 FASTQ
        r2_path: Path to R2 FASTQ
        output_dir: Output directory
        threads: Number of threads

    Returns:
        Dict with 'success', 'merged_path', 'unmerged_r1', 'unmerged_r2'
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    merged_path = output_dir / "merged.fastq.gz"
    unmerged_r1 = output_dir / "unmerged_R1.fastq.gz"
    unmerged_r2 = output_dir / "unmerged_R2.fastq.gz"

    # Use BBMerge
    cmd = [
        "bbmerge.sh",
        f"in1={r1_path}",
        f"in2={r2_path}",
        f"out={merged_path}",
        f"outu1={unmerged_r1}",
        f"outu2={unmerged_r2}",
        f"threads={threads}",
        "mininsert=35",
        "minoverlap=15",
        "strict=t",
    ]

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=600,
        )

        if result.returncode == 0 and merged_path.exists():
            return {
                'success': True,
                'merged_path': merged_path,
                'unmerged_r1': unmerged_r1 if unmerged_r1.exists() else None,
                'unmerged_r2': unmerged_r2 if unmerged_r2.exists() else None,
            }
        else:
            return {
                'success': False,
                'error': result.stderr or "BBMerge failed",
            }

    except FileNotFoundError:
        logger.warning("BBMerge not found, skipping merge step")
        return {
            'success': False,
            'error': "BBMerge not installed",
        }
    except subprocess.TimeoutExpired:
        return {
            'success': False,
            'error': "BBMerge timed out",
        }
    except Exception as e:
        return {
            'success': False,
            'error': str(e),
        }


def run_collapse(
    input_path: Path,
    output_dir: Path,
) -> Dict:
    """
    Collapse identical sequences, keeping one representative.

    Args:
        input_path: Path to input FASTQ
        output_dir: Output directory

    Returns:
        Dict with 'success', 'collapsed_path', 'unique_count', 'total_count'
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    collapsed_path = output_dir / "collapsed.fastq.gz"

    # Track unique sequences: seq -> (header, qual, count)
    unique_seqs: Dict[str, Tuple[str, str, int]] = {}
    total = 0

    open_func = gzip.open if str(input_path).endswith('.gz') else open

    try:
        with open_func(input_path, 'rt') as f:
            while True:
                header = f.readline().strip()
                if not header:
                    break
                seq = f.readline().strip()
                plus = f.readline().strip()
                qual = f.readline().strip()

                total += 1

                if seq in unique_seqs:
                    # Keep higher quality version
                    old_header, old_qual, count = unique_seqs[seq]
                    old_score = sum(ord(c) - 33 for c in old_qual)
                    new_score = sum(ord(c) - 33 for c in qual)

                    if new_score > old_score:
                        unique_seqs[seq] = (header, qual, count + 1)
                    else:
                        unique_seqs[seq] = (old_header, old_qual, count + 1)
                else:
                    unique_seqs[seq] = (header, qual, 1)

        # Write collapsed file
        with gzip.open(collapsed_path, 'wt') as out:
            for seq, (header, qual, count) in unique_seqs.items():
                # Append count to header for reference
                new_header = f"{header};count={count}"
                out.write(f"{new_header}\n{seq}\n+\n{qual}\n")

        return {
            'success': True,
            'collapsed_path': collapsed_path,
            'total_count': total,
            'unique_count': len(unique_seqs),
        }

    except Exception as e:
        logger.error(f"Collapse failed: {e}")
        return {
            'success': False,
            'error': str(e),
        }


def count_reads(path: Path) -> int:
    """Count reads in a FASTQ file."""
    if not path or not path.exists():
        return 0

    count = 0
    open_func = gzip.open if str(path).endswith('.gz') else open

    try:
        with open_func(path, 'rt') as f:
            for line in f:
                if line.startswith('@'):
                    count += 1
    except Exception:
        pass

    return count
