"""
K-mer based HDR/WT classification.

Uses short k-mers (typically 12-mers) spanning the edit site to rapidly
classify reads without alignment.

Author: Kevin R. Roy
"""

import gzip
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Dict, List, Optional, Set

from ..utils.sequence import (
    extract_unique_kmers,
    generate_edit_kmers,
    reverse_complement,
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

        For barcode-style templates (same flanking sequences, different middle),
        generates k-mers that span the barcode boundaries to discriminate between
        templates.

        Args:
            reference: Wild-type reference sequence
            hdr_templates: Dict mapping template_id → HDR template sequence
            edit_positions_by_template: Dict mapping template_id → list of edit positions
            contamination_sequence: Optional sequence for contamination k-mers
            kmer_size: Size of k-mers to use

        Returns:
            MultiTemplateKmerClassifier instance
        """
        # Detect if this is a barcode-style template set
        is_barcode_style = cls._detect_barcode_templates(hdr_templates)

        if is_barcode_style:
            return cls._from_barcode_templates(
                reference, hdr_templates, kmer_size, contamination_sequence
            )
        else:
            return cls._from_general_templates(
                reference, hdr_templates, edit_positions_by_template,
                kmer_size, contamination_sequence
            )

    @classmethod
    def _detect_barcode_templates(cls, hdr_templates: Dict[str, str]) -> bool:
        """Check if templates have barcode-style structure (shared flanks)."""
        if len(hdr_templates) < 2:
            return False

        templates = list(hdr_templates.values())
        first = templates[0].upper()

        # Find common prefix length
        min_prefix = len(first)
        for t in templates[1:]:
            t_upper = t.upper()
            common = 0
            for i in range(min(len(first), len(t_upper))):
                if first[i] == t_upper[i]:
                    common += 1
                else:
                    break
            min_prefix = min(min_prefix, common)

        # Find common suffix length
        min_suffix = len(first)
        for t in templates[1:]:
            t_upper = t.upper()
            common = 0
            for i in range(1, min(len(first), len(t_upper)) + 1):
                if first[-i] == t_upper[-i]:
                    common += 1
                else:
                    break
            min_suffix = min(min_suffix, common)

        # Barcode style if significant common prefix and suffix
        return min_prefix >= 20 and min_suffix >= 20

    @classmethod
    def _from_barcode_templates(
        cls,
        reference: str,
        hdr_templates: Dict[str, str],
        kmer_size: int,
        contamination_sequence: Optional[str],
    ) -> 'MultiTemplateKmerClassifier':
        """
        Create classifier for barcode-style templates.

        For barcoded experiments, k-mers must SPAN the entire barcode plus
        flanking bases to ensure uniqueness. K-mer size is automatically
        increased to barcode_size + buffer (minimum 10bp buffer).

        Example: For 10bp barcodes with 20bp k-mers (10bp + 10bp buffer):
        - 9bp_left + 10bp_barcode + 1bp_right
        - 5bp_left + 10bp_barcode + 5bp_right
        - 1bp_left + 10bp_barcode + 9bp_right
        """
        templates = list(hdr_templates.values())
        first = templates[0].upper()

        # Find common prefix and suffix lengths (shared homology arms)
        prefix_len = len(first)
        suffix_len = len(first)

        for t in templates[1:]:
            t_upper = t.upper()
            # Update prefix
            common = 0
            for i in range(min(prefix_len, len(t_upper))):
                if first[i] == t_upper[i]:
                    common += 1
                else:
                    break
            prefix_len = common

            # Update suffix
            common = 0
            for i in range(1, min(suffix_len, len(t_upper)) + 1):
                if first[-i] == t_upper[-i]:
                    common += 1
                else:
                    break
            suffix_len = common

        # Barcode region in template coordinates
        barcode_start = prefix_len
        barcode_len = len(first) - prefix_len - suffix_len
        barcode_end = barcode_start + barcode_len

        # For barcode templates, k-mer size must span the entire barcode + buffer
        # Minimum k-mer size = barcode_len + 10 (5bp on each side minimum)
        min_kmer_size = barcode_len + 10
        effective_kmer_size = max(kmer_size, min_kmer_size)

        # Find where templates align in reference
        first_anchor = first[:min(30, prefix_len)]
        ref_upper = reference.upper()
        ref_offset = ref_upper.find(first_anchor)
        if ref_offset == -1:
            ref_offset = 0

        # WT k-mers: k-mers from reference spanning where barcode would be inserted
        all_wt_kmers: Set[str] = set()
        ref_barcode_pos = ref_offset + barcode_start

        for start in range(max(0, ref_barcode_pos - effective_kmer_size + 1),
                          min(len(reference) - effective_kmer_size + 1, ref_barcode_pos + 1)):
            kmer = reference[start:start + effective_kmer_size].upper()
            all_wt_kmers.add(kmer)
            all_wt_kmers.add(reverse_complement(kmer))

        # HDR k-mers per template: k-mers that SPAN the entire barcode
        hdr_kmers_by_template: Dict[str, Set[str]] = {}

        for template_id, template_seq in hdr_templates.items():
            template_upper = template_seq.upper()
            hdr_kmers: Set[str] = set()

            # Generate k-mers that fully span the barcode
            # K-mer must contain entire barcode: start <= barcode_start AND end >= barcode_end
            # With k-mer size K and barcode length B:
            # - Earliest start: barcode_end - K (so k-mer ends at barcode_end)
            # - Latest start: barcode_start (so k-mer starts at barcode_start)
            # Valid starts: max(0, barcode_end - K) to min(barcode_start, len - K)

            earliest_start = max(0, barcode_end - effective_kmer_size)
            latest_start = min(barcode_start, len(template_upper) - effective_kmer_size)

            for start in range(earliest_start, latest_start + 1):
                end = start + effective_kmer_size
                # Verify this k-mer spans the entire barcode
                if start <= barcode_start and end >= barcode_end:
                    kmer = template_upper[start:end]
                    hdr_kmers.add(kmer)
                    hdr_kmers.add(reverse_complement(kmer))

            # Remove any WT k-mers (shouldn't be any, but be safe)
            hdr_kmers_by_template[template_id] = hdr_kmers - all_wt_kmers

        # Verify uniqueness - if k-mers are shared, increase k-mer size dynamically
        # Continue increasing until uniqueness is achieved or max size is reached
        all_hdr = []
        for kmers in hdr_kmers_by_template.values():
            all_hdr.extend(kmers)
        from collections import Counter
        kmer_counts = Counter(all_hdr)
        shared_count = sum(1 for c in kmer_counts.values() if c > 1)

        # Maximum k-mer size: template length minus small buffer for flanking
        min_template_len = min(len(t) for t in hdr_templates.values())
        max_kmer_size = min_template_len - 4  # Leave 2bp on each side

        if shared_count > 0 and effective_kmer_size < max_kmer_size:
            # Dynamically increase k-mer size until uniqueness achieved
            import logging
            logger = logging.getLogger(__name__)
            logger.debug(
                f"K-mer size {effective_kmer_size}bp has {shared_count} shared k-mers, "
                f"increasing to {effective_kmer_size + 2}bp"
            )
            return cls._from_barcode_templates(
                reference, hdr_templates, effective_kmer_size + 2, contamination_sequence
            )

        # Warn if uniqueness still not achieved at max size
        if shared_count > 0:
            import logging
            logger = logging.getLogger(__name__)
            logger.warning(
                f"Could not achieve k-mer uniqueness even at {effective_kmer_size}bp. "
                f"{shared_count} k-mers are shared between templates. "
                f"Classification may have reduced accuracy."
            )

        # Contamination k-mers
        contamination_kmers: Set[str] = set()
        if contamination_sequence:
            contamination_kmers = extract_unique_kmers(
                contamination_sequence,
                kmer_size=effective_kmer_size,
                exclude_sequences=[reference] + list(hdr_templates.values()),
            )

        return cls(
            wt_kmers=all_wt_kmers,
            hdr_kmers_by_template=hdr_kmers_by_template,
            contamination_kmers=contamination_kmers,
            kmer_size=effective_kmer_size,
        )

    @classmethod
    def _from_general_templates(
        cls,
        reference: str,
        hdr_templates: Dict[str, str],
        edit_positions_by_template: Dict[str, List[int]],
        kmer_size: int,
        contamination_sequence: Optional[str],
    ) -> 'MultiTemplateKmerClassifier':
        """Create classifier using general edit position approach."""
        all_wt_kmers: Set[str] = set()
        hdr_kmers_by_template: Dict[str, Set[str]] = {}

        for template_id, template_seq in hdr_templates.items():
            edit_positions = edit_positions_by_template.get(template_id, [])

            if not edit_positions:
                hdr_kmers_by_template[template_id] = set()
                continue

            wt_kmers, hdr_kmers = generate_edit_kmers(
                reference, template_seq, edit_positions, kmer_size
            )
            all_wt_kmers.update(wt_kmers)
            hdr_only = hdr_kmers - all_wt_kmers
            hdr_kmers_by_template[template_id] = hdr_only

        contamination_kmers: Set[str] = set()
        if contamination_sequence:
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
