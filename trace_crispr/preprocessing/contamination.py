"""
Contamination filtering using k-mer based detection.

Author: Kevin R. Roy
"""

import gzip
from pathlib import Path
from typing import Set, Optional, Tuple
from dataclasses import dataclass
import logging

from ..utils.sequence import extract_unique_kmers, reverse_complement

logger = logging.getLogger(__name__)


@dataclass
class ContaminationFilterResult:
    """Result from contamination filtering."""
    success: bool
    output_r1: Path
    output_r2: Optional[Path]
    total_reads: int = 0
    filtered_reads: int = 0
    contamination_rate: float = 0.0
    error_message: Optional[str] = None


def filter_contamination_fastq(
    r1_fastq: Path,
    r2_fastq: Optional[Path],
    output_r1: Path,
    output_r2: Optional[Path],
    contamination_kmers: Set[str],
) -> ContaminationFilterResult:
    """
    Filter reads containing contamination k-mers.

    For paired-end reads, if EITHER R1 or R2 contains a contamination
    k-mer, the entire pair is removed.

    Args:
        r1_fastq: Input R1 FASTQ
        r2_fastq: Optional input R2 FASTQ
        output_r1: Output R1 FASTQ path
        output_r2: Optional output R2 FASTQ path
        contamination_kmers: Set of contamination k-mers

    Returns:
        ContaminationFilterResult with statistics
    """
    if not contamination_kmers:
        # No filtering needed
        import shutil
        shutil.copy(r1_fastq, output_r1)
        if r2_fastq and output_r2:
            shutil.copy(r2_fastq, output_r2)

        return ContaminationFilterResult(
            success=True,
            output_r1=output_r1,
            output_r2=output_r2,
        )

    try:
        total_reads = 0
        filtered_reads = 0

        open_func = gzip.open if str(r1_fastq).endswith('.gz') else open
        out_func = gzip.open if str(output_r1).endswith('.gz') else open

        if r2_fastq:
            # Paired-end filtering
            with open_func(r1_fastq, 'rt') as f1, \
                 open_func(r2_fastq, 'rt') as f2, \
                 out_func(output_r1, 'wt') as o1, \
                 out_func(output_r2, 'wt') as o2:

                while True:
                    # Read R1 record
                    header1 = f1.readline()
                    if not header1:
                        break
                    seq1 = f1.readline()
                    plus1 = f1.readline()
                    qual1 = f1.readline()

                    # Read R2 record
                    header2 = f2.readline()
                    seq2 = f2.readline()
                    plus2 = f2.readline()
                    qual2 = f2.readline()

                    total_reads += 1

                    # Check for contamination in either read
                    has_contamination = (
                        _contains_kmer(seq1.strip().upper(), contamination_kmers) or
                        _contains_kmer(seq2.strip().upper(), contamination_kmers)
                    )

                    if has_contamination:
                        filtered_reads += 1
                    else:
                        # Write passing reads
                        o1.write(header1 + seq1 + plus1 + qual1)
                        o2.write(header2 + seq2 + plus2 + qual2)
        else:
            # Single-end filtering
            with open_func(r1_fastq, 'rt') as f1, \
                 out_func(output_r1, 'wt') as o1:

                while True:
                    header = f1.readline()
                    if not header:
                        break
                    seq = f1.readline()
                    plus = f1.readline()
                    qual = f1.readline()

                    total_reads += 1

                    if _contains_kmer(seq.strip().upper(), contamination_kmers):
                        filtered_reads += 1
                    else:
                        o1.write(header + seq + plus + qual)

        contamination_rate = filtered_reads / total_reads if total_reads > 0 else 0

        return ContaminationFilterResult(
            success=True,
            output_r1=output_r1,
            output_r2=output_r2,
            total_reads=total_reads,
            filtered_reads=filtered_reads,
            contamination_rate=contamination_rate,
        )

    except Exception as e:
        return ContaminationFilterResult(
            success=False,
            output_r1=output_r1,
            output_r2=output_r2,
            error_message=str(e),
        )


def _contains_kmer(sequence: str, kmers: Set[str]) -> bool:
    """Check if sequence contains any k-mer from the set."""
    if not kmers:
        return False

    # Get k-mer size from first k-mer
    k = len(next(iter(kmers)))

    for i in range(len(sequence) - k + 1):
        if sequence[i:i+k] in kmers:
            return True

    return False


def create_contamination_filter(
    contamination_fasta: Path,
    reference_fasta: Optional[Path] = None,
    kmer_size: int = 12,
) -> Set[str]:
    """
    Create a set of contamination k-mers for filtering.

    Args:
        contamination_fasta: Path to contamination sequence FASTA
        reference_fasta: Optional reference FASTA to exclude k-mers from
        kmer_size: K-mer size

    Returns:
        Set of contamination k-mers (including reverse complements)
    """
    # Read contamination sequence
    contamination_seq = _read_fasta(contamination_fasta)

    # Read reference sequence if provided
    exclude_seqs = []
    if reference_fasta:
        exclude_seqs.append(_read_fasta(reference_fasta))

    # Generate unique k-mers
    kmers = extract_unique_kmers(
        contamination_seq,
        kmer_size=kmer_size,
        exclude_sequences=exclude_seqs,
    )

    logger.info(f"Created {len(kmers)} contamination k-mers")

    return kmers


def _read_fasta(path: Path) -> str:
    """Read first sequence from FASTA file."""
    sequence = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if sequence:
                    break
                continue
            sequence.append(line.upper())
    return ''.join(sequence)
