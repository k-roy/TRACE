"""
K-mer based HDR/WT classification.

Uses short k-mers (typically 12-mers) spanning the edit site to rapidly
classify reads without alignment.

Author: Kevin R. Roy
"""

from dataclasses import dataclass
from typing import Set, Tuple, Optional, List
from enum import Enum
import gzip
from pathlib import Path

from ..utils.sequence import (
    reverse_complement,
    generate_edit_kmers,
    extract_unique_kmers,
)


class KmerClassification(Enum):
    """K-mer based classification categories."""
    WILD_TYPE = 'wt'
    HDR = 'hdr'
    CONTAMINATION = 'contamination'
    AMBIGUOUS = 'ambiguous'
    UNKNOWN = 'unknown'


@dataclass
class KmerClassifier:
    """K-mer based read classifier."""
    wt_kmers: Set[str]
    hdr_kmers: Set[str]
    contamination_kmers: Set[str]
    kmer_size: int

    @classmethod
    def from_sequences(
        cls,
        reference: str,
        hdr_template: str,
        edit_positions: List[int],
        contamination_sequence: Optional[str] = None,
        kmer_size: int = 12,
    ) -> 'KmerClassifier':
        """
        Create a KmerClassifier from sequences.

        Args:
            reference: Wild-type reference sequence
            hdr_template: HDR template sequence (same length as reference)
            edit_positions: List of positions where edits occur
            contamination_sequence: Optional sequence for contamination k-mers
            kmer_size: Size of k-mers to use

        Returns:
            KmerClassifier instance
        """
        # Generate WT and HDR k-mers spanning edit sites
        wt_kmers, hdr_kmers = generate_edit_kmers(
            reference, hdr_template, edit_positions, kmer_size
        )

        # Generate contamination k-mers if provided
        contamination_kmers = set()
        if contamination_sequence:
            # Get k-mers unique to contamination sequence
            contamination_kmers = extract_unique_kmers(
                contamination_sequence,
                kmer_size=kmer_size,
                exclude_sequences=[reference, hdr_template],
            )

        return cls(
            wt_kmers=wt_kmers,
            hdr_kmers=hdr_kmers,
            contamination_kmers=contamination_kmers,
            kmer_size=kmer_size,
        )

    def classify_sequence(self, sequence: str) -> KmerClassification:
        """
        Classify a sequence based on k-mer content.

        Args:
            sequence: DNA sequence to classify

        Returns:
            KmerClassification enum value
        """
        seq_upper = sequence.upper()

        # Check for matches in each category
        has_wt = self._contains_any_kmer(seq_upper, self.wt_kmers)
        has_hdr = self._contains_any_kmer(seq_upper, self.hdr_kmers)
        has_contamination = self._contains_any_kmer(seq_upper, self.contamination_kmers)

        # Contamination takes priority
        if has_contamination and not has_hdr:
            return KmerClassification.CONTAMINATION

        # Check for unambiguous classification
        if has_hdr and not has_wt:
            return KmerClassification.HDR
        if has_wt and not has_hdr:
            return KmerClassification.WILD_TYPE
        if has_wt and has_hdr:
            return KmerClassification.AMBIGUOUS

        return KmerClassification.UNKNOWN

    def classify_read_pair(
        self,
        r1_sequence: str,
        r2_sequence: str
    ) -> KmerClassification:
        """
        Classify a paired-end read based on k-mer content.

        If either read contains HDR k-mers, the pair is classified as HDR.
        This increases sensitivity for HDR detection.

        Args:
            r1_sequence: R1 DNA sequence
            r2_sequence: R2 DNA sequence

        Returns:
            KmerClassification enum value
        """
        r1_class = self.classify_sequence(r1_sequence)
        r2_class = self.classify_sequence(r2_sequence)

        # If either is contamination, check if the other is HDR
        if r1_class == KmerClassification.CONTAMINATION:
            if r2_class == KmerClassification.HDR:
                return KmerClassification.HDR
            return KmerClassification.CONTAMINATION
        if r2_class == KmerClassification.CONTAMINATION:
            if r1_class == KmerClassification.HDR:
                return KmerClassification.HDR
            return KmerClassification.CONTAMINATION

        # HDR takes priority (if either read has HDR, classify as HDR)
        if r1_class == KmerClassification.HDR or r2_class == KmerClassification.HDR:
            return KmerClassification.HDR

        # WT if both are WT
        if r1_class == KmerClassification.WILD_TYPE and r2_class == KmerClassification.WILD_TYPE:
            return KmerClassification.WILD_TYPE

        # If one is WT and other is unknown, classify as WT
        if r1_class == KmerClassification.WILD_TYPE or r2_class == KmerClassification.WILD_TYPE:
            return KmerClassification.WILD_TYPE

        # Ambiguous if any is ambiguous
        if r1_class == KmerClassification.AMBIGUOUS or r2_class == KmerClassification.AMBIGUOUS:
            return KmerClassification.AMBIGUOUS

        return KmerClassification.UNKNOWN

    def _contains_any_kmer(self, sequence: str, kmers: Set[str]) -> bool:
        """Check if sequence contains any k-mer from the set."""
        if not kmers:
            return False

        for i in range(len(sequence) - self.kmer_size + 1):
            kmer = sequence[i:i + self.kmer_size]
            if kmer in kmers:
                return True

        return False


@dataclass
class KmerResults:
    """Results from k-mer classification of a sample."""
    total_reads: int
    wt_count: int
    hdr_count: int
    contamination_count: int
    ambiguous_count: int
    unknown_count: int

    @property
    def wt_rate(self) -> float:
        return self.wt_count / self.total_reads if self.total_reads > 0 else 0

    @property
    def hdr_rate(self) -> float:
        return self.hdr_count / self.total_reads if self.total_reads > 0 else 0

    @property
    def contamination_rate(self) -> float:
        return self.contamination_count / self.total_reads if self.total_reads > 0 else 0

    @property
    def classifiable_count(self) -> int:
        return self.wt_count + self.hdr_count

    @property
    def classifiable_hdr_rate(self) -> float:
        """HDR rate among classifiable (WT + HDR) reads."""
        total = self.classifiable_count
        return self.hdr_count / total if total > 0 else 0


def classify_fastq_kmer(
    r1_path: Path,
    r2_path: Optional[Path],
    classifier: KmerClassifier,
) -> KmerResults:
    """
    Classify all reads in a FASTQ file using k-mer approach.

    Args:
        r1_path: Path to R1 FASTQ file (gzipped or not)
        r2_path: Optional path to R2 FASTQ file
        classifier: KmerClassifier instance

    Returns:
        KmerResults with classification counts
    """
    counts = {
        KmerClassification.WILD_TYPE: 0,
        KmerClassification.HDR: 0,
        KmerClassification.CONTAMINATION: 0,
        KmerClassification.AMBIGUOUS: 0,
        KmerClassification.UNKNOWN: 0,
    }

    # Open files
    open_func = gzip.open if str(r1_path).endswith('.gz') else open

    if r2_path:
        # Paired-end mode
        with open_func(r1_path, 'rt') as f1, open_func(r2_path, 'rt') as f2:
            while True:
                # Read R1 record
                header1 = f1.readline().strip()
                if not header1:
                    break
                seq1 = f1.readline().strip()
                f1.readline()  # +
                f1.readline()  # qual

                # Read R2 record
                f2.readline()  # header
                seq2 = f2.readline().strip()
                f2.readline()  # +
                f2.readline()  # qual

                classification = classifier.classify_read_pair(seq1, seq2)
                counts[classification] += 1
    else:
        # Single-end mode
        with open_func(r1_path, 'rt') as f1:
            while True:
                header = f1.readline().strip()
                if not header:
                    break
                seq = f1.readline().strip()
                f1.readline()  # +
                f1.readline()  # qual

                classification = classifier.classify_sequence(seq)
                counts[classification] += 1

    total = sum(counts.values())

    return KmerResults(
        total_reads=total,
        wt_count=counts[KmerClassification.WILD_TYPE],
        hdr_count=counts[KmerClassification.HDR],
        contamination_count=counts[KmerClassification.CONTAMINATION],
        ambiguous_count=counts[KmerClassification.AMBIGUOUS],
        unknown_count=counts[KmerClassification.UNKNOWN],
    )
