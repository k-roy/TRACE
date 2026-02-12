"""
Adapter trimming using BBDuk.

Author: Kevin R. Roy
"""

import subprocess
from pathlib import Path
from typing import Optional, Tuple
from dataclasses import dataclass
import logging

logger = logging.getLogger(__name__)

# Standard adapter sequences
TN5_ADAPTERS = [
    "CTGTCTCTTATACACATCT",
    "AGATGTGTATAAGAGACAG",  # Reverse complement
]

TRUSEQ_ADAPTERS = [
    "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",  # Read 1
    "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",  # Read 2
]


@dataclass
class TrimmingResult:
    """Result from adapter trimming."""
    success: bool
    output_r1: Path
    output_r2: Optional[Path]
    reads_in: int = 0
    reads_out: int = 0
    bases_trimmed: int = 0
    error_message: Optional[str] = None


def run_bbduk_trim(
    r1_fastq: Path,
    r2_fastq: Optional[Path],
    output_r1: Path,
    output_r2: Optional[Path],
    adapters: list,
    min_length: int = 50,
    ktrim: str = 'r',
    k: int = 23,
    mink: int = 11,
    hdist: int = 1,
    threads: int = 4,
) -> TrimmingResult:
    """
    Run BBDuk for adapter trimming.

    Args:
        r1_fastq: Input R1 FASTQ
        r2_fastq: Optional input R2 FASTQ
        output_r1: Output R1 FASTQ path
        output_r2: Optional output R2 FASTQ path
        adapters: List of adapter sequences
        min_length: Minimum read length after trimming
        ktrim: Trim mode ('r' = right, 'l' = left)
        k: Kmer length for matching
        mink: Minimum kmer length at ends
        hdist: Hamming distance for matching
        threads: Number of threads

    Returns:
        TrimmingResult with trimming statistics
    """
    try:
        # Write adapters to temporary file
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
            for i, adapter in enumerate(adapters):
                f.write(f">adapter_{i}\n{adapter}\n")
            adapter_file = Path(f.name)

        # Build command
        cmd = [
            'bbduk.sh',
            f'in={r1_fastq}',
            f'out={output_r1}',
            f'ref={adapter_file}',
            f'ktrim={ktrim}',
            f'k={k}',
            f'mink={mink}',
            f'hdist={hdist}',
            f'minlen={min_length}',
            f'threads={threads}',
            'tpe=t',  # Trim pairs evenly
            'tbo=t',  # Trim by overlap
        ]

        if r2_fastq and output_r2:
            cmd.extend([
                f'in2={r2_fastq}',
                f'out2={output_r2}',
            ])

        # Run BBDuk
        result = subprocess.run(cmd, capture_output=True, check=True)

        # Parse stats from stderr
        stats = result.stderr.decode()
        reads_in = 0
        reads_out = 0

        for line in stats.split('\n'):
            if 'Input:' in line:
                parts = line.split()
                for i, p in enumerate(parts):
                    if p == 'reads':
                        try:
                            reads_in = int(parts[i-1])
                        except (ValueError, IndexError):
                            pass
            if 'Result:' in line:
                parts = line.split()
                for i, p in enumerate(parts):
                    if p == 'reads':
                        try:
                            reads_out = int(parts[i-1])
                        except (ValueError, IndexError):
                            pass

        # Clean up adapter file
        adapter_file.unlink()

        return TrimmingResult(
            success=True,
            output_r1=output_r1,
            output_r2=output_r2,
            reads_in=reads_in,
            reads_out=reads_out,
        )

    except subprocess.CalledProcessError as e:
        return TrimmingResult(
            success=False,
            output_r1=output_r1,
            output_r2=output_r2,
            error_message=e.stderr.decode() if e.stderr else str(e),
        )
    except Exception as e:
        return TrimmingResult(
            success=False,
            output_r1=output_r1,
            output_r2=output_r2,
            error_message=str(e),
        )


def trim_adapters(
    r1_fastq: Path,
    r2_fastq: Optional[Path],
    output_dir: Path,
    library_type: str = 'auto',
    threads: int = 4,
) -> TrimmingResult:
    """
    Trim adapters based on library type.

    Args:
        r1_fastq: Input R1 FASTQ
        r2_fastq: Optional input R2 FASTQ
        output_dir: Output directory
        library_type: 'Tn5', 'TruSeq', 'both', or 'auto'
        threads: Number of threads

    Returns:
        TrimmingResult
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Determine adapters
    if library_type == 'Tn5':
        adapters = TN5_ADAPTERS
    elif library_type == 'TruSeq':
        adapters = TRUSEQ_ADAPTERS
    elif library_type in ('both', 'auto'):
        adapters = TN5_ADAPTERS + TRUSEQ_ADAPTERS
    else:
        adapters = TN5_ADAPTERS + TRUSEQ_ADAPTERS

    # Create output paths
    r1_name = r1_fastq.name.replace('.fastq', '_trimmed.fastq').replace('.fq', '_trimmed.fq')
    output_r1 = output_dir / r1_name

    output_r2 = None
    if r2_fastq:
        r2_name = r2_fastq.name.replace('.fastq', '_trimmed.fastq').replace('.fq', '_trimmed.fq')
        output_r2 = output_dir / r2_name

    logger.info(f"Trimming adapters (library type: {library_type})")

    return run_bbduk_trim(
        r1_fastq=r1_fastq,
        r2_fastq=r2_fastq,
        output_r1=output_r1,
        output_r2=output_r2,
        adapters=adapters,
        threads=threads,
    )
