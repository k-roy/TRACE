"""
Editing outcome classification for CRISPR experiments.

Classifies reads into:
- HDR (Perfect): Exact match to HDR template
- HDR (Imperfect): HDR signature with additional edits
- NHEJ Indel: Indel at cut site without HDR signature
- Large Deletion: Deletion >= threshold overlapping target region
- Wild-type: Matches reference exactly
- Unclassified: None of the above

Author: Kevin R. Roy
"""

from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass, field
from enum import Enum
from collections import Counter
import pysam

from .cigar import (
    get_deletions_from_cigar,
    get_insertions_from_cigar,
    count_mismatches_in_region,
)


class EditingOutcome(Enum):
    """Editing outcome categories."""
    HDR_PERFECT = 'hdr_perfect'
    HDR_IMPERFECT = 'hdr_imperfect'
    NHEJ_INSERTION = 'nhej_insertion'
    NHEJ_DELETION = 'nhej_deletion'
    LARGE_DELETION = 'large_deletion'
    WILD_TYPE = 'wild_type'
    UNCLASSIFIED = 'unclassified'
    UNMAPPED = 'unmapped'


@dataclass
class ClassificationResult:
    """Result of read classification."""
    outcome: EditingOutcome
    details: Dict = field(default_factory=dict)
    confidence: float = 1.0
    qc_flags: List[str] = field(default_factory=list)


def get_hdr_signature_positions(wt_seq: str, hdr_seq: str) -> List[Tuple[int, str, str]]:
    """
    Identify positions where HDR differs from WT.

    Args:
        wt_seq: Wild-type sequence
        hdr_seq: HDR template sequence

    Returns:
        List of (position, wt_base, hdr_base) tuples
    """
    if len(wt_seq) != len(hdr_seq):
        raise ValueError("WT and HDR sequences must be same length")

    signature_positions = []
    for i, (wt_base, hdr_base) in enumerate(zip(wt_seq, hdr_seq)):
        if wt_base.upper() != hdr_base.upper():
            signature_positions.append((i, wt_base.upper(), hdr_base.upper()))

    return signature_positions


def count_hdr_signature_matches(
    read: pysam.AlignedSegment,
    hdr_signature: List[Tuple[int, str, str]],
    ref_offset: int = 0
) -> Tuple[int, int, List[int]]:
    """
    Count how many HDR signature positions match in a read.

    Args:
        read: pysam AlignedSegment object
        hdr_signature: List of (position, wt_base, hdr_base) from get_hdr_signature_positions
        ref_offset: Offset to convert signature positions to reference coordinates

    Returns:
        Tuple of (n_hdr_matches, n_positions_covered, matched_positions)
    """
    if read.query_sequence is None:
        return 0, 0, []

    # Get aligned pairs
    aligned_pairs = read.get_aligned_pairs(with_seq=True)

    # Build lookup of ref_pos -> query_base
    ref_to_query = {}
    for query_pos, ref_pos, ref_base in aligned_pairs:
        if ref_pos is not None and query_pos is not None:
            ref_to_query[ref_pos] = read.query_sequence[query_pos].upper()

    n_matches = 0
    n_covered = 0
    matched_positions = []

    for sig_pos, wt_base, hdr_base in hdr_signature:
        ref_pos = sig_pos + ref_offset

        if ref_pos in ref_to_query:
            n_covered += 1
            query_base = ref_to_query[ref_pos]

            if query_base == hdr_base:
                n_matches += 1
                matched_positions.append(sig_pos)

    return n_matches, n_covered, matched_positions


def check_indels_at_cut_site(
    read: pysam.AlignedSegment,
    cut_site: int,
    window: int = 5
) -> Tuple[List[Dict], List[Dict]]:
    """
    Check for indels near the cut site.

    Args:
        read: pysam AlignedSegment object
        cut_site: Position of CRISPR cut site in reference
        window: Window around cut site to check

    Returns:
        Tuple of (deletions, insertions) found near cut site
    """
    target_region = (cut_site - window, cut_site + window)

    deletions = get_deletions_from_cigar(read, min_size=1, target_region=target_region)
    insertions = get_insertions_from_cigar(read, min_size=1, target_region=target_region)

    del_dicts = [{'start': d.ref_start, 'size': d.size} for d in deletions]
    ins_dicts = [{'position': i.ref_position, 'size': i.size, 'seq': i.inserted_seq}
                 for i in insertions]

    return del_dicts, ins_dicts


def is_wild_type_match(
    read: pysam.AlignedSegment,
    cut_site: int,
    window_size: int = 20,
    max_mismatches: int = 0
) -> bool:
    """
    Check if read matches wild-type reference around cut site.

    Args:
        read: pysam AlignedSegment object
        cut_site: Position of CRISPR cut site
        window_size: Window around cut site to check
        max_mismatches: Maximum allowed mismatches

    Returns:
        True if read matches WT within tolerance
    """
    region_start = cut_site - window_size
    region_end = cut_site + window_size

    # Check for any indels in region
    deletions = get_deletions_from_cigar(read, min_size=1,
                                         target_region=(region_start, region_end))
    insertions = get_insertions_from_cigar(read, min_size=1,
                                           target_region=(region_start, region_end))

    if deletions or insertions:
        return False

    # Check mismatches
    n_matches, n_mismatches, _ = count_mismatches_in_region(
        read, region_start, region_end
    )

    return n_mismatches <= max_mismatches


def classify_read(
    read: pysam.AlignedSegment,
    hdr_signature: List[Tuple[int, str, str]],
    cut_site: int,
    ref_offset: int = 0,
    window_size: int = 20,
    large_del_min_size: int = 21,
    hdr_threshold: float = 0.8,
) -> ClassificationResult:
    """
    Classify a read into an editing outcome category.

    Classification hierarchy:
    1. Unmapped -> UNMAPPED
    2. Large deletion (>= large_del_min_size bp) -> LARGE_DELETION
    3. HDR signature present without indels -> HDR_PERFECT
    4. HDR signature present with indels -> HDR_IMPERFECT
    5. Insertion at cut site -> NHEJ_INSERTION
    6. Deletion at cut site -> NHEJ_DELETION
    7. Matches WT -> WILD_TYPE
    8. Otherwise -> UNCLASSIFIED

    Args:
        read: pysam AlignedSegment object
        hdr_signature: HDR signature positions from get_hdr_signature_positions
        cut_site: Position of CRISPR cut site in reference
        ref_offset: Offset to convert signature positions to reference coordinates
        window_size: Window around cut site for analysis
        large_del_min_size: Minimum size for large deletion
        hdr_threshold: Fraction of HDR positions that must match

    Returns:
        ClassificationResult with outcome and details
    """
    qc_flags = []
    details = {}

    # Check if unmapped
    if read.is_unmapped:
        return ClassificationResult(
            outcome=EditingOutcome.UNMAPPED,
            details={},
            confidence=1.0,
            qc_flags=[]
        )

    # Define target region for large deletion detection
    target_region = (cut_site - window_size - 50, cut_site + window_size + 50)

    # Check for large deletions first
    large_dels = get_deletions_from_cigar(
        read,
        min_size=large_del_min_size,
        target_region=target_region
    )

    if large_dels:
        # Get largest deletion
        largest_del = max(large_dels, key=lambda d: d.size)

        details['deletions'] = [{'start': d.ref_start, 'size': d.size} for d in large_dels]
        details['largest_deletion_size'] = largest_del.size

        return ClassificationResult(
            outcome=EditingOutcome.LARGE_DELETION,
            details=details,
            confidence=0.9,
            qc_flags=qc_flags
        )

    # Check HDR signature
    n_hdr_matches, n_covered, matched_positions = count_hdr_signature_matches(
        read, hdr_signature, ref_offset
    )

    hdr_match_fraction = n_hdr_matches / len(hdr_signature) if hdr_signature else 0
    details['hdr_match_fraction'] = hdr_match_fraction
    details['hdr_matched_positions'] = matched_positions

    # Check for indels at cut site
    cut_site_deletions, cut_site_insertions = check_indels_at_cut_site(
        read, cut_site, window=5
    )

    has_deletions = bool(cut_site_deletions)
    has_insertions = bool(cut_site_insertions)
    details['cut_site_deletions'] = cut_site_deletions
    details['cut_site_insertions'] = cut_site_insertions

    # Classification logic
    if hdr_match_fraction >= hdr_threshold:
        if not has_deletions and not has_insertions:
            return ClassificationResult(
                outcome=EditingOutcome.HDR_PERFECT,
                details=details,
                confidence=hdr_match_fraction,
                qc_flags=qc_flags
            )
        else:
            return ClassificationResult(
                outcome=EditingOutcome.HDR_IMPERFECT,
                details=details,
                confidence=hdr_match_fraction * 0.9,
                qc_flags=qc_flags
            )

    if has_insertions:
        return ClassificationResult(
            outcome=EditingOutcome.NHEJ_INSERTION,
            details=details,
            confidence=0.9,
            qc_flags=qc_flags
        )

    if has_deletions:
        return ClassificationResult(
            outcome=EditingOutcome.NHEJ_DELETION,
            details=details,
            confidence=0.9,
            qc_flags=qc_flags
        )

    # Check if wild-type
    if is_wild_type_match(read, cut_site, window_size):
        return ClassificationResult(
            outcome=EditingOutcome.WILD_TYPE,
            details=details,
            confidence=0.95,
            qc_flags=qc_flags
        )

    return ClassificationResult(
        outcome=EditingOutcome.UNCLASSIFIED,
        details=details,
        confidence=0.5,
        qc_flags=['no_clear_classification']
    )


def summarize_classifications(
    classifications: List[ClassificationResult]
) -> Dict:
    """
    Summarize a list of classifications into counts and rates.

    Args:
        classifications: List of ClassificationResult objects

    Returns:
        Dict with counts and rates for each outcome
    """
    counts = Counter(c.outcome for c in classifications)

    total = len(classifications)
    summary = {
        'total_reads': total,
    }

    for outcome in EditingOutcome:
        count = counts.get(outcome, 0)
        summary[f'{outcome.value}_count'] = count
        summary[f'{outcome.value}_rate'] = count / total if total > 0 else 0

    # Combined HDR rate
    hdr_total = counts.get(EditingOutcome.HDR_PERFECT, 0) + counts.get(EditingOutcome.HDR_IMPERFECT, 0)
    summary['hdr_total_count'] = hdr_total
    summary['hdr_total_rate'] = hdr_total / total if total > 0 else 0

    # Combined NHEJ rate
    nhej_total = (
        counts.get(EditingOutcome.NHEJ_INSERTION, 0) +
        counts.get(EditingOutcome.NHEJ_DELETION, 0) +
        counts.get(EditingOutcome.LARGE_DELETION, 0)
    )
    summary['nhej_total_count'] = nhej_total
    summary['nhej_total_rate'] = nhej_total / total if total > 0 else 0

    # Editing rate (anything not WT or unmapped)
    non_edited = counts.get(EditingOutcome.WILD_TYPE, 0) + counts.get(EditingOutcome.UNMAPPED, 0)
    edited = total - non_edited
    summary['edited_count'] = edited
    summary['edited_rate'] = edited / total if total > 0 else 0

    return summary
