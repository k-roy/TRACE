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

from collections import Counter
from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Tuple

import pysam

from .cigar import (
    count_mismatches_in_region,
    get_deletions_from_cigar,
    get_insertions_from_cigar,
)


class EditingOutcome(Enum):
    """
    Comprehensive editing outcome categories.

    HDR categories:
    - HDR_COMPLETE: All donor-encoded edits present, no other modifications
    - HDR_PARTIAL: ≥1 but not all donor-encoded edits present
    - HDR_PLUS_NHEJ_INDEL: HDR edits + classical NHEJ indel (0-2bp microhomology)
    - HDR_PLUS_MMEJ_INDEL: HDR edits + MMEJ indel (>2bp microhomology)
    - HDR_PLUS_OTHER: HDR edits + modifications NOT at cut site

    DONOR_CAPTURE: HDR edits present + extra donor sequence duplicated at site
                   (mutually exclusive with HDR_* - takes precedence)

    Non-HDR repair outcomes:
    - NHEJ_INDEL: Classical NHEJ indel at cut site (0-2bp microhomology)
    - MMEJ_INDEL: Microhomology-mediated indel at cut site (>2bp microhomology)

    Other:
    - WT: Wild-type / unedited
    - NON_DONOR_SNV: SNVs not matching donor, no indels
    - UNCLASSIFIED: Does not fit other categories
    - UNMAPPED: Read did not align
    """
    # HDR categories
    HDR_COMPLETE = 'hdr_complete'
    HDR_PARTIAL = 'hdr_partial'
    HDR_PLUS_NHEJ_INDEL = 'hdr_plus_nhej_indel'
    HDR_PLUS_MMEJ_INDEL = 'hdr_plus_mmej_indel'
    HDR_PLUS_OTHER = 'hdr_plus_other'

    # Donor capture (mutually exclusive with HDR_*)
    DONOR_CAPTURE = 'donor_capture'

    # Non-HDR repair outcomes
    NHEJ_INDEL = 'nhej_indel'
    MMEJ_INDEL = 'mmej_indel'

    # Other outcomes
    WT = 'wt'
    NON_DONOR_SNV = 'non_donor_snv'
    UNCLASSIFIED = 'unclassified'
    UNMAPPED = 'unmapped'


# Microhomology threshold for NHEJ vs MMEJ classification
# ≤2bp = classical NHEJ, >2bp = MMEJ (alt-NHEJ)
MMEJ_MICROHOMOLOGY_THRESHOLD = 2


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


def calculate_deletion_microhomology(
    ref_seq: str,
    deletion_start: int,
    deletion_end: int,
    max_search: int = 20
) -> Tuple[int, int, str, str]:
    """
    Calculate microhomology at deletion boundaries.

    Microhomology indicates the deletion repair mechanism - NHEJ-mediated
    deletions often show microhomology at the breakpoint junctions.

    Args:
        ref_seq: Reference sequence
        deletion_start: Start position of deletion in reference
        deletion_end: End position of deletion in reference
        max_search: Maximum distance to search for microhomology (default: 20bp)

    Returns:
        Tuple of (left_mh_length, right_mh_length, left_mh_seq, right_mh_seq)
        where mh = microhomology
    """
    ref_upper = ref_seq.upper()
    deletion_length = deletion_end - deletion_start
    deleted_seq = ref_upper[deletion_start:deletion_end] if deletion_end <= len(ref_upper) else ""

    if not deleted_seq:
        return 0, 0, "", ""

    # How far left could this deletion be shifted? (left microhomology)
    left_shift = 0
    for i in range(1, min(max_search, deletion_length, deletion_start) + 1):
        # Check if prefix of deleted seq matches suffix before deletion
        if deletion_start - i >= 0 and deleted_seq[:i] == ref_upper[deletion_start - i:deletion_start]:
            left_shift = i
        # NO BREAK - continue to find maximum shift (handles interrupted repeats)

    # How far right could this deletion be shifted? (right microhomology)
    right_shift = 0
    for i in range(1, min(max_search, deletion_length, len(ref_upper) - deletion_end) + 1):
        # Check if suffix of deleted seq matches prefix after deletion
        if deletion_end + i <= len(ref_upper) and deleted_seq[-i:] == ref_upper[deletion_end:deletion_end + i]:
            right_shift = i
        # NO BREAK - continue to find maximum shift

    # Extract the actual microhomology sequences
    left_mh_seq = ref_upper[deletion_start - left_shift:deletion_start] if left_shift > 0 else ""
    right_mh_seq = ref_upper[deletion_end:deletion_end + right_shift] if right_shift > 0 else ""

    return left_shift, right_shift, left_mh_seq, right_mh_seq


def classify_indel_mechanism(
    ref_seq: str,
    deletion_start: int,
    deletion_end: int,
    mh_threshold: int = MMEJ_MICROHOMOLOGY_THRESHOLD
) -> Tuple[str, int]:
    """
    Classify indel repair mechanism based on microhomology.

    Args:
        ref_seq: Reference sequence
        deletion_start: Start position of deletion
        deletion_end: End position of deletion
        mh_threshold: Microhomology threshold (default: 2bp)
                     ≤threshold = NHEJ, >threshold = MMEJ

    Returns:
        Tuple of (mechanism: 'nhej' or 'mmej', max_microhomology_length)
    """
    if not ref_seq:
        # Without reference sequence, default to NHEJ
        return 'nhej', 0

    left_mh, right_mh, _, _ = calculate_deletion_microhomology(
        ref_seq, deletion_start, deletion_end
    )
    max_mh = max(left_mh, right_mh)

    if max_mh > mh_threshold:
        return 'mmej', max_mh
    else:
        return 'nhej', max_mh


def detect_donor_capture(
    read: pysam.AlignedSegment,
    donor_seq: str,
    cut_site: int,
    min_match_length: int = 30,
    window: int = 100
) -> Tuple[bool, Dict]:
    """
    Detect donor capture: extra donor sequence duplicated at target site.

    Donor capture occurs when HDR happens but additional donor sequence
    is inserted beyond the intended edits.

    Args:
        read: pysam AlignedSegment object
        donor_seq: Full donor template sequence
        cut_site: Position of CRISPR cut site
        min_match_length: Minimum length of donor match to call capture
        window: Window around cut site to search

    Returns:
        Tuple of (is_donor_capture, details_dict)
    """
    if not donor_seq or read.query_sequence is None:
        return False, {}

    # Get insertions near cut site
    target_region = (cut_site - window, cut_site + window)
    insertions = get_insertions_from_cigar(read, min_size=min_match_length, target_region=target_region)

    if not insertions:
        return False, {}

    donor_upper = donor_seq.upper()

    for ins in insertions:
        ins_seq = ins.inserted_seq.upper()
        if len(ins_seq) < min_match_length:
            continue

        # Check if insertion matches donor sequence
        if ins_seq in donor_upper:
            return True, {
                'insertion_position': ins.ref_position,
                'insertion_size': ins.size,
                'insertion_seq': ins.inserted_seq,
                'donor_match': True
            }

        # Also check reverse complement for antisense donors
        rc_donor = _reverse_complement(donor_upper)
        if ins_seq in rc_donor:
            return True, {
                'insertion_position': ins.ref_position,
                'insertion_size': ins.size,
                'insertion_seq': ins.inserted_seq,
                'donor_match': True,
                'antisense': True
            }

    return False, {}


def _reverse_complement(seq: str) -> str:
    """Return reverse complement of DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq.upper()))


def check_non_donor_snvs(
    read: pysam.AlignedSegment,
    hdr_signature: List[Tuple[int, str, str]],
    cut_site: int,
    window_size: int = 20,
    ref_offset: int = 0
) -> Tuple[bool, List[Dict]]:
    """
    Check for SNVs that don't match the donor template.

    Args:
        read: pysam AlignedSegment
        hdr_signature: HDR signature positions
        cut_site: Cut site position
        window_size: Window around cut site
        ref_offset: Reference offset

    Returns:
        Tuple of (has_non_donor_snvs, list of SNV details)
    """
    if read.query_sequence is None:
        return False, []

    region_start = cut_site - window_size
    region_end = cut_site + window_size

    # Get all mismatches in region
    aligned_pairs = read.get_aligned_pairs(with_seq=True)

    # Build set of HDR signature positions
    hdr_positions = {pos + ref_offset for pos, _, _ in hdr_signature}

    non_donor_snvs = []
    for query_pos, ref_pos, ref_base in aligned_pairs:
        if ref_pos is None or query_pos is None:
            continue
        if ref_pos < region_start or ref_pos > region_end:
            continue

        query_base = read.query_sequence[query_pos].upper()
        if ref_base and ref_base.upper() != query_base:
            # This is a mismatch - check if it's an HDR signature position
            if ref_pos not in hdr_positions:
                non_donor_snvs.append({
                    'position': ref_pos,
                    'ref_base': ref_base.upper() if ref_base else 'N',
                    'query_base': query_base
                })

    return len(non_donor_snvs) > 0, non_donor_snvs


def classify_read(
    read: pysam.AlignedSegment,
    hdr_signature: List[Tuple[int, str, str]],
    cut_site: int,
    ref_seq: str = "",
    donor_seq: str = "",
    ref_offset: int = 0,
    window_size: int = 20,
    large_del_min_size: int = 50,
    hdr_threshold: float = 0.8,
    partial_hdr_threshold: float = 0.0,  # Any HDR edit = partial
) -> ClassificationResult:
    """
    Classify a read into an editing outcome category.

    Classification hierarchy:
    1. UNMAPPED: Read did not align
    2. Check for HDR signature (donor edits):
       a. If HDR edits + extra donor sequence → DONOR_CAPTURE
       b. If ALL HDR edits + no indels → HDR_COMPLETE
       c. If ALL HDR edits + indel → HDR_PLUS_NHEJ_INDEL or HDR_PLUS_MMEJ_INDEL
       d. If SOME HDR edits → HDR_PARTIAL
       e. If HDR edits + non-cut-site edits → HDR_PLUS_OTHER
    3. No HDR edits:
       a. Indel at cut site with >2bp microhomology → MMEJ_INDEL
       b. Indel at cut site with ≤2bp microhomology → NHEJ_INDEL
       c. SNVs not matching donor → NON_DONOR_SNV
       d. No edits → WT
       e. Otherwise → UNCLASSIFIED

    Args:
        read: pysam AlignedSegment object
        hdr_signature: HDR signature positions from get_hdr_signature_positions
        cut_site: Position of CRISPR cut site in reference
        ref_seq: Reference sequence (for microhomology analysis)
        donor_seq: Full donor template sequence (for donor capture detection)
        ref_offset: Offset to convert signature positions to reference coordinates
        window_size: Window around cut site for analysis
        large_del_min_size: Minimum size for large deletion flag
        hdr_threshold: Fraction of HDR positions for HDR_COMPLETE (default: 0.8)
        partial_hdr_threshold: Fraction for HDR_PARTIAL (default: 0.0 = any edit)

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

    # Define target region
    target_region = (cut_site - window_size - 50, cut_site + window_size + 50)

    # ========== Gather all edit information ==========

    # Check HDR signature
    n_hdr_matches, n_covered, matched_positions = count_hdr_signature_matches(
        read, hdr_signature, ref_offset
    )

    n_hdr_total = len(hdr_signature) if hdr_signature else 0
    hdr_match_fraction = n_hdr_matches / n_hdr_total if n_hdr_total > 0 else 0
    has_any_hdr_edits = n_hdr_matches > 0
    has_all_hdr_edits = hdr_match_fraction >= hdr_threshold

    details['hdr_match_fraction'] = hdr_match_fraction
    details['hdr_matched_positions'] = matched_positions
    details['hdr_n_matches'] = n_hdr_matches
    details['hdr_n_total'] = n_hdr_total

    # Check for indels at cut site (±5bp window)
    cut_site_deletions, cut_site_insertions = check_indels_at_cut_site(
        read, cut_site, window=5
    )

    has_cut_site_deletions = bool(cut_site_deletions)
    has_cut_site_insertions = bool(cut_site_insertions)
    has_cut_site_indels = has_cut_site_deletions or has_cut_site_insertions

    details['cut_site_deletions'] = cut_site_deletions
    details['cut_site_insertions'] = cut_site_insertions

    # Check for large deletions (for metadata flag)
    large_dels = get_deletions_from_cigar(read, min_size=large_del_min_size, target_region=target_region)
    cut_spanning_large_dels = [d for d in large_dels if d.ref_start <= cut_site <= d.ref_end]
    is_large_deletion = bool(cut_spanning_large_dels)
    details['is_large_deletion'] = is_large_deletion

    if is_large_deletion:
        largest_del = max(cut_spanning_large_dels, key=lambda d: d.size)
        details['largest_deletion_size'] = largest_del.size
        details['deletions'] = [{'start': d.ref_start, 'end': d.ref_end, 'size': d.size}
                                for d in cut_spanning_large_dels]

    # Determine indel repair mechanism if there are indels
    indel_mechanism = None
    max_microhomology = 0

    if has_cut_site_deletions and ref_seq:
        # Use the largest deletion for mechanism classification
        all_dels = get_deletions_from_cigar(read, min_size=1, target_region=(cut_site - 5, cut_site + 5))
        if all_dels:
            largest_cut_del = max(all_dels, key=lambda d: d.size)
            indel_mechanism, max_microhomology = classify_indel_mechanism(
                ref_seq, largest_cut_del.ref_start, largest_cut_del.ref_end
            )
            details['indel_mechanism'] = indel_mechanism
            details['max_microhomology'] = max_microhomology

            # Add microhomology details
            left_mh, right_mh, left_seq, right_seq = calculate_deletion_microhomology(
                ref_seq, largest_cut_del.ref_start, largest_cut_del.ref_end
            )
            details['microhomology_left'] = left_mh
            details['microhomology_right'] = right_mh
            details['microhomology_left_seq'] = left_seq
            details['microhomology_right_seq'] = right_seq

    elif has_cut_site_insertions:
        # Insertions are classified as NHEJ (blunt insertions)
        indel_mechanism = 'nhej'
        details['indel_mechanism'] = 'nhej'

    # Check for donor capture (HDR + extra donor sequence)
    is_donor_capture = False
    donor_capture_details = {}
    if has_any_hdr_edits and donor_seq:
        is_donor_capture, donor_capture_details = detect_donor_capture(
            read, donor_seq, cut_site
        )
        if is_donor_capture:
            details['donor_capture'] = donor_capture_details

    # Check for non-donor SNVs
    has_non_donor_snvs, non_donor_snv_list = check_non_donor_snvs(
        read, hdr_signature, cut_site, window_size, ref_offset
    )
    if has_non_donor_snvs:
        details['non_donor_snvs'] = non_donor_snv_list

    # ========== Classification Logic ==========

    # 1. DONOR_CAPTURE takes precedence (mutually exclusive with HDR_*)
    if is_donor_capture:
        return ClassificationResult(
            outcome=EditingOutcome.DONOR_CAPTURE,
            details=details,
            confidence=0.85,
            qc_flags=qc_flags
        )

    # 2. HDR categories
    if has_any_hdr_edits:
        if has_all_hdr_edits:
            # All HDR edits present
            if not has_cut_site_indels:
                # HDR_COMPLETE: All edits, no additional indels
                return ClassificationResult(
                    outcome=EditingOutcome.HDR_COMPLETE,
                    details=details,
                    confidence=hdr_match_fraction,
                    qc_flags=qc_flags
                )
            else:
                # HDR + indel at cut site
                if indel_mechanism == 'mmej':
                    return ClassificationResult(
                        outcome=EditingOutcome.HDR_PLUS_MMEJ_INDEL,
                        details=details,
                        confidence=hdr_match_fraction * 0.9,
                        qc_flags=qc_flags
                    )
                else:
                    return ClassificationResult(
                        outcome=EditingOutcome.HDR_PLUS_NHEJ_INDEL,
                        details=details,
                        confidence=hdr_match_fraction * 0.9,
                        qc_flags=qc_flags
                    )
        else:
            # Partial HDR (some but not all edits)
            if has_non_donor_snvs and not has_cut_site_indels:
                # HDR edits + non-donor modifications not at cut site
                return ClassificationResult(
                    outcome=EditingOutcome.HDR_PLUS_OTHER,
                    details=details,
                    confidence=hdr_match_fraction * 0.8,
                    qc_flags=qc_flags
                )
            else:
                return ClassificationResult(
                    outcome=EditingOutcome.HDR_PARTIAL,
                    details=details,
                    confidence=hdr_match_fraction * 0.9,
                    qc_flags=qc_flags
                )

    # 3. No HDR edits - check for non-HDR repair outcomes
    if has_cut_site_indels:
        if indel_mechanism == 'mmej':
            return ClassificationResult(
                outcome=EditingOutcome.MMEJ_INDEL,
                details=details,
                confidence=0.9,
                qc_flags=qc_flags
            )
        else:
            return ClassificationResult(
                outcome=EditingOutcome.NHEJ_INDEL,
                details=details,
                confidence=0.9,
                qc_flags=qc_flags
            )

    # 4. Check for non-donor SNVs (no indels, no HDR)
    if has_non_donor_snvs:
        return ClassificationResult(
            outcome=EditingOutcome.NON_DONOR_SNV,
            details=details,
            confidence=0.7,
            qc_flags=['snv_only']
        )

    # 5. Check if wild-type
    if is_wild_type_match(read, cut_site, window_size):
        return ClassificationResult(
            outcome=EditingOutcome.WT,
            details=details,
            confidence=0.95,
            qc_flags=qc_flags
        )

    # 6. Unclassified
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
        Dict with counts and rates for each outcome, including aggregated metrics:
        - hdr_total: All HDR categories combined
        - nhej_mmej_total: NHEJ_INDEL + MMEJ_INDEL
        - edited_total: All edited reads (not WT or unmapped)
    """
    counts = Counter(c.outcome for c in classifications)

    total = len(classifications)
    summary = {
        'total_reads': total,
    }

    # Per-category counts and rates
    for outcome in EditingOutcome:
        count = counts.get(outcome, 0)
        summary[f'{outcome.value}_count'] = count
        summary[f'{outcome.value}_rate'] = count / total if total > 0 else 0

    # ========== Aggregated HDR metrics ==========
    # HDR_total = HDR_COMPLETE + HDR_PARTIAL + HDR_PLUS_NHEJ + HDR_PLUS_MMEJ + HDR_PLUS_OTHER
    hdr_total = (
        counts.get(EditingOutcome.HDR_COMPLETE, 0) +
        counts.get(EditingOutcome.HDR_PARTIAL, 0) +
        counts.get(EditingOutcome.HDR_PLUS_NHEJ_INDEL, 0) +
        counts.get(EditingOutcome.HDR_PLUS_MMEJ_INDEL, 0) +
        counts.get(EditingOutcome.HDR_PLUS_OTHER, 0)
    )
    summary['hdr_total_count'] = hdr_total
    summary['hdr_total_rate'] = hdr_total / total if total > 0 else 0

    # ========== Aggregated NHEJ/MMEJ metrics ==========
    # NHEJ_MMEJ_total = NHEJ_INDEL + MMEJ_INDEL
    nhej_mmej_total = (
        counts.get(EditingOutcome.NHEJ_INDEL, 0) +
        counts.get(EditingOutcome.MMEJ_INDEL, 0)
    )
    summary['nhej_mmej_total_count'] = nhej_mmej_total
    summary['nhej_mmej_total_rate'] = nhej_mmej_total / total if total > 0 else 0

    # ========== Edited total ==========
    # Edited_total = HDR_total + NHEJ_MMEJ_total + DONOR_CAPTURE
    # (excludes WT, NON_DONOR_SNV, UNCLASSIFIED, UNMAPPED)
    edited_total = hdr_total + nhej_mmej_total + counts.get(EditingOutcome.DONOR_CAPTURE, 0)
    summary['edited_total_count'] = edited_total
    summary['edited_total_rate'] = edited_total / total if total > 0 else 0

    # ========== Large deletion flag summary ==========
    # Count reads with is_large_deletion flag
    n_large_deletions = sum(
        1 for c in classifications
        if c.details.get('is_large_deletion', False)
    )
    summary['large_deletion_flag_count'] = n_large_deletions
    summary['large_deletion_flag_rate'] = n_large_deletions / total if total > 0 else 0

    # ========== Mechanism breakdown for indels ==========
    n_mmej_mechanism = sum(
        1 for c in classifications
        if c.details.get('indel_mechanism') == 'mmej'
    )
    n_nhej_mechanism = sum(
        1 for c in classifications
        if c.details.get('indel_mechanism') == 'nhej'
    )
    summary['mmej_mechanism_count'] = n_mmej_mechanism
    summary['nhej_mechanism_count'] = n_nhej_mechanism

    return summary
