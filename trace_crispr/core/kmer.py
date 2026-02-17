"""
K-mer based HDR/WT classification.

Uses short k-mers (typically 12-mers) spanning the edit site to rapidly
classify reads without alignment.

Author: Kevin R. Roy
"""

from dataclasses import dataclass, field
from typing import Set, Tuple, Optional, List, Dict
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


@dataclass
class MultiTemplateKmerResult:
    """K-mer classification result for multi-template analysis."""
    is_wt: bool
    is_contamination: bool
    template_matches: Dict[str, int]  # template_id → k-mer match count
    best_template: Optional[str]  # template with most matches (None if WT/contam)
    is_ambiguous: bool  # True if multiple templates match equally

    @property
    def has_hdr(self) -> bool:
        """True if any template has matches."""
        return any(count > 0 for count in self.template_matches.values())

    @property
    def classification(self) -> str:
        """Get classification as string for output."""
        if self.is_contamination:
            return "CONTAMINATION"
        if self.is_wt and not self.has_hdr:
            return "WT"
        if self.is_ambiguous:
            return "AMBIGUOUS"
        if self.best_template:
            return f"HDR:{self.best_template}"
        return "UNKNOWN"


@dataclass
class MultiTemplateKmerClassifier:
    """
    K-mer classifier that distinguishes between multiple HDR templates.

    Used for barcode screening experiments where we want to detect which
    barcode (if any) is present in each read, and identify unexpected
    barcodes for purity checking.
    """
    wt_kmers: Set[str]
    hdr_kmers_by_template: Dict[str, Set[str]]  # template_id → kmers
    contamination_kmers: Set[str]
    kmer_size: int

    @classmethod
    def from_multi_template(
        cls,
        reference: str,
        hdr_templates: Dict[str, str],
        edit_positions_by_template: Dict[str, List[int]],
        contamination_sequence: Optional[str] = None,
        kmer_size: int = 12,
    ) -> 'MultiTemplateKmerClassifier':
        """
        Create classifier from multiple HDR templates.

        Args:
            reference: Wild-type reference sequence
            hdr_templates: Dict mapping template_id → HDR template sequence
            edit_positions_by_template: Dict mapping template_id → list of edit positions
            contamination_sequence: Optional sequence for contamination k-mers
            kmer_size: Size of k-mers to use

        Returns:
            MultiTemplateKmerClassifier instance
        """
        # Generate WT k-mers (shared across templates)
        # Use edit positions from first template to get representative WT k-mers
        all_wt_kmers: Set[str] = set()
        hdr_kmers_by_template: Dict[str, Set[str]] = {}

        for template_id, template_seq in hdr_templates.items():
            edit_positions = edit_positions_by_template.get(template_id, [])

            if not edit_positions:
                # No edits detected, skip this template
                hdr_kmers_by_template[template_id] = set()
                continue

            # Generate WT and HDR k-mers for this template
            wt_kmers, hdr_kmers = generate_edit_kmers(
                reference, template_seq, edit_positions, kmer_size
            )

            # Add WT k-mers to shared set
            all_wt_kmers.update(wt_kmers)

            # Remove any WT k-mers that might appear in this template's HDR set
            # (shouldn't happen normally, but be safe)
            hdr_only = hdr_kmers - all_wt_kmers
            hdr_kmers_by_template[template_id] = hdr_only

        # Remove any HDR k-mers that appear in multiple templates
        # These would cause ambiguous classifications
        # (Actually, we want to keep them - the classifier will handle ambiguity)

        # Generate contamination k-mers if provided
        contamination_kmers: Set[str] = set()
        if contamination_sequence:
            # Get k-mers unique to contamination
            all_hdr_kmers = set()
            for kmers in hdr_kmers_by_template.values():
                all_hdr_kmers.update(kmers)

            contamination_kmers = extract_unique_kmers(
                contamination_sequence,
                kmer_size=kmer_size,
                exclude_sequences=[reference] + list(hdr_templates.values()),
            )

        return cls(
            wt_kmers=all_wt_kmers,
            hdr_kmers_by_template=hdr_kmers_by_template,
            contamination_kmers=contamination_kmers,
            kmer_size=kmer_size,
        )

    def classify_sequence(self, sequence: str) -> MultiTemplateKmerResult:
        """
        Classify a sequence against all templates.

        Returns:
            MultiTemplateKmerResult with:
            - is_wt: True if sequence has WT k-mers but no HDR k-mers
            - is_contamination: True if contamination k-mers detected
            - template_matches: Dict of template_id → match count
            - best_template: Template with most matches (or None)
            - is_ambiguous: True if multiple templates tie for most matches
        """
        seq_upper = sequence.upper()

        # Check for contamination
        has_contamination = self._contains_any_kmer(seq_upper, self.contamination_kmers)

        # Check for WT k-mers
        has_wt = self._contains_any_kmer(seq_upper, self.wt_kmers)

        # Count matches for each template
        template_matches: Dict[str, int] = {}
        for template_id, hdr_kmers in self.hdr_kmers_by_template.items():
            count = self._count_kmer_matches(seq_upper, hdr_kmers)
            template_matches[template_id] = count

        # Determine best template
        max_count = max(template_matches.values()) if template_matches else 0
        best_templates = [
            tid for tid, count in template_matches.items()
            if count == max_count and count > 0
        ]

        # Handle ambiguity
        is_ambiguous = len(best_templates) > 1
        best_template = best_templates[0] if len(best_templates) == 1 else None

        # Contamination takes priority if no HDR detected
        if has_contamination and max_count == 0:
            return MultiTemplateKmerResult(
                is_wt=False,
                is_contamination=True,
                template_matches=template_matches,
                best_template=None,
                is_ambiguous=False,
            )

        # Determine if this is WT (has WT k-mers but no HDR)
        is_wt = has_wt and max_count == 0

        return MultiTemplateKmerResult(
            is_wt=is_wt,
            is_contamination=has_contamination and not is_wt and max_count == 0,
            template_matches=template_matches,
            best_template=best_template,
            is_ambiguous=is_ambiguous,
        )

    def classify_read_pair(
        self,
        r1_sequence: str,
        r2_sequence: str
    ) -> MultiTemplateKmerResult:
        """
        Classify a paired-end read based on k-mer content.

        Combines results from both reads, summing template matches.
        """
        r1_result = self.classify_sequence(r1_sequence)
        r2_result = self.classify_sequence(r2_sequence)

        # Combine template matches
        combined_matches: Dict[str, int] = {}
        all_templates = set(r1_result.template_matches.keys()) | set(r2_result.template_matches.keys())
        for tid in all_templates:
            combined_matches[tid] = (
                r1_result.template_matches.get(tid, 0) +
                r2_result.template_matches.get(tid, 0)
            )

        # Determine best template from combined
        max_count = max(combined_matches.values()) if combined_matches else 0
        best_templates = [
            tid for tid, count in combined_matches.items()
            if count == max_count and count > 0
        ]

        is_ambiguous = len(best_templates) > 1
        best_template = best_templates[0] if len(best_templates) == 1 else None

        # Combine WT and contamination status
        has_wt = r1_result.is_wt or r2_result.is_wt
        has_contamination = r1_result.is_contamination or r2_result.is_contamination

        # If either has HDR, classify as HDR
        is_wt = has_wt and max_count == 0

        return MultiTemplateKmerResult(
            is_wt=is_wt,
            is_contamination=has_contamination and not is_wt and max_count == 0,
            template_matches=combined_matches,
            best_template=best_template,
            is_ambiguous=is_ambiguous,
        )

    def _contains_any_kmer(self, sequence: str, kmers: Set[str]) -> bool:
        """Check if sequence contains any k-mer from the set."""
        if not kmers:
            return False
        for i in range(len(sequence) - self.kmer_size + 1):
            kmer = sequence[i:i + self.kmer_size]
            if kmer in kmers:
                return True
        return False

    def _count_kmer_matches(self, sequence: str, kmers: Set[str]) -> int:
        """Count how many k-mers from the set are found in the sequence."""
        if not kmers:
            return 0
        count = 0
        for i in range(len(sequence) - self.kmer_size + 1):
            kmer = sequence[i:i + self.kmer_size]
            if kmer in kmers:
                count += 1
        return count


@dataclass
class MultiTemplateKmerResults:
    """Results from multi-template k-mer classification of a sample."""
    total_reads: int
    wt_count: int
    hdr_counts_by_template: Dict[str, int]
    contamination_count: int
    ambiguous_count: int
    unknown_count: int

    @property
    def total_hdr_count(self) -> int:
        return sum(self.hdr_counts_by_template.values())

    @property
    def wt_rate(self) -> float:
        return self.wt_count / self.total_reads if self.total_reads > 0 else 0

    @property
    def total_hdr_rate(self) -> float:
        return self.total_hdr_count / self.total_reads if self.total_reads > 0 else 0

    def hdr_rate_for_template(self, template_id: str) -> float:
        count = self.hdr_counts_by_template.get(template_id, 0)
        return count / self.total_reads if self.total_reads > 0 else 0


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


def classify_fastq_multi_template(
    r1_path: Path,
    r2_path: Optional[Path],
    classifier: MultiTemplateKmerClassifier,
) -> MultiTemplateKmerResults:
    """
    Classify all reads in a FASTQ file using multi-template k-mer approach.

    Args:
        r1_path: Path to R1 FASTQ file (gzipped or not)
        r2_path: Optional path to R2 FASTQ file
        classifier: MultiTemplateKmerClassifier instance

    Returns:
        MultiTemplateKmerResults with per-template classification counts
    """
    # Use a dict for mutable count tracking
    counts = {
        'wt': 0,
        'contamination': 0,
        'ambiguous': 0,
        'unknown': 0,
    }
    hdr_counts: Dict[str, int] = {tid: 0 for tid in classifier.hdr_kmers_by_template.keys()}

    def update_counts(result: MultiTemplateKmerResult):
        """Update counts based on classification result."""
        if result.is_wt:
            counts['wt'] += 1
        elif result.is_contamination:
            counts['contamination'] += 1
        elif result.is_ambiguous:
            counts['ambiguous'] += 1
        elif result.best_template:
            hdr_counts[result.best_template] = hdr_counts.get(result.best_template, 0) + 1
        else:
            counts['unknown'] += 1

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

                result = classifier.classify_read_pair(seq1, seq2)
                update_counts(result)
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

                result = classifier.classify_sequence(seq)
                update_counts(result)

    total = counts['wt'] + counts['contamination'] + counts['ambiguous'] + counts['unknown'] + sum(hdr_counts.values())

    return MultiTemplateKmerResults(
        total_reads=total,
        wt_count=counts['wt'],
        hdr_counts_by_template=hdr_counts,
        contamination_count=counts['contamination'],
        ambiguous_count=counts['ambiguous'],
        unknown_count=counts['unknown'],
    )
