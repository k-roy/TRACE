"""
Edit-distance based HDR detection.

For each mismatch/indel in a read aligned to WT reference:
- Check if incorporating it makes the read closer to the donor
- If yes, classify as donor-encoded edit
- If no, classify as error or NHEJ

This approach:
- Handles any number of SNVs without combinatorial explosion
- Works naturally with insertions and deletions
- Doesn't require pre-computing all variant combinations

Author: Kevin R. Roy
"""

import logging
from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)


def needleman_wunsch(seq1: str, seq2: str, match=2, mismatch=-1, gap=-2) -> Tuple[str, str, int]:
    """
    Global sequence alignment using Needleman-Wunsch algorithm.

    Returns aligned sequences with gaps represented as '-'.

    Args:
        seq1, seq2: Sequences to align
        match: Score for matching bases
        mismatch: Penalty for mismatches
        gap: Penalty for gaps (insertions/deletions)

    Returns:
        (aligned_seq1, aligned_seq2, score)

    Raises:
        ValueError: If either sequence is empty
    """
    if not seq1 or not seq2:
        raise ValueError(
            f"Cannot align empty sequences. seq1 length={len(seq1) if seq1 else 0}, "
            f"seq2 length={len(seq2) if seq2 else 0}"
        )
    n, m = len(seq1), len(seq2)

    # Initialize scoring matrix
    score = [[0] * (m + 1) for _ in range(n + 1)]

    # Initialize first row and column (gap penalties)
    for i in range(n + 1):
        score[i][0] = gap * i
    for j in range(m + 1):
        score[0][j] = gap * j

    # Fill scoring matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match_score = score[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
            delete_score = score[i-1][j] + gap
            insert_score = score[i][j-1] + gap
            score[i][j] = max(match_score, delete_score, insert_score)

    # Traceback to get alignment
    align1, align2 = [], []
    i, j = n, m

    while i > 0 or j > 0:
        if i > 0 and j > 0:
            current_score = score[i][j]
            diag_score = score[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)

            if current_score == diag_score:
                align1.append(seq1[i-1])
                align2.append(seq2[j-1])
                i -= 1
                j -= 1
            elif current_score == score[i-1][j] + gap:
                align1.append(seq1[i-1])
                align2.append('-')
                i -= 1
            else:
                align1.append('-')
                align2.append(seq2[j-1])
                j -= 1
        elif i > 0:
            align1.append(seq1[i-1])
            align2.append('-')
            i -= 1
        else:
            align1.append('-')
            align2.append(seq2[j-1])
            j -= 1

    return ''.join(reversed(align1)), ''.join(reversed(align2)), score[n][m]


def smith_waterman(seq1: str, seq2: str, match=2, mismatch=-1, gap=-2) -> Tuple[str, str, int, int, int]:
    """
    Local sequence alignment using Smith-Waterman algorithm.

    Finds the best local alignment between seq1 and seq2, ignoring terminal mismatches.
    This is ideal for aligning donor templates with 5'/3' overhangs to reference sequences.

    Args:
        seq1: Reference sequence
        seq2: Query sequence (e.g., donor template)
        match: Score for matching bases
        mismatch: Penalty for mismatches
        gap: Penalty for gaps

    Returns:
        (aligned_seq1, aligned_seq2, score, start_pos_in_seq1, start_pos_in_seq2)

    Raises:
        ValueError: If either sequence is empty
    """
    if not seq1 or not seq2:
        raise ValueError(
            f"Cannot align empty sequences. seq1 length={len(seq1) if seq1 else 0}, "
            f"seq2 length={len(seq2) if seq2 else 0}"
        )
    n, m = len(seq1), len(seq2)

    # Initialize scoring matrix (local alignment allows 0)
    score = [[0] * (m + 1) for _ in range(n + 1)]

    # Track maximum score and position
    max_score = 0
    max_i, max_j = 0, 0

    # Fill scoring matrix (first row/column stay 0 for local alignment)
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match_score = score[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
            delete_score = score[i-1][j] + gap
            insert_score = score[i][j-1] + gap
            score[i][j] = max(0, match_score, delete_score, insert_score)

            if score[i][j] > max_score:
                max_score = score[i][j]
                max_i, max_j = i, j

    # Traceback from maximum score position
    align1, align2 = [], []
    i, j = max_i, max_j
    start_i, start_j = i, j

    while i > 0 and j > 0 and score[i][j] > 0:
        current_score = score[i][j]
        diag_score = score[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)

        if current_score == diag_score:
            align1.append(seq1[i-1])
            align2.append(seq2[j-1])
            start_i, start_j = i-1, j-1
            i -= 1
            j -= 1
        elif current_score == score[i-1][j] + gap:
            align1.append(seq1[i-1])
            align2.append('-')
            start_i = i-1
            i -= 1
        elif current_score == score[i][j-1] + gap:
            align1.append('-')
            align2.append(seq2[j-1])
            start_j = j-1
            j -= 1
        else:
            break

    return ''.join(reversed(align1)), ''.join(reversed(align2)), max_score, start_i, start_j


class EditType(Enum):
    """Type of edit detected in a read."""
    SNV = "snv"
    INSERTION = "insertion"
    DELETION = "deletion"


@dataclass
class DetectedEdit:
    """A single edit detected in a read."""
    edit_type: EditType
    ref_position: int  # Position in reference
    ref_base: str  # Base(s) in reference (empty for insertion)
    read_base: str  # Base(s) in read (empty for deletion)
    is_donor_encoded: bool  # Does this edit move toward donor?
    distance_to_cut: int = 0

    @property
    def size(self) -> int:
        """Size of the edit (1 for SNV, length for indel)."""
        if self.edit_type == EditType.SNV:
            return 1
        elif self.edit_type == EditType.INSERTION:
            return len(self.read_base)
        else:  # DELETION
            return len(self.ref_base)


@dataclass
class DonorSignature:
    """Complete signature of donor-encoded edits including SNVs and indels."""
    snvs: Dict[int, Tuple[str, str]]  # position -> (ref_base, donor_base)
    insertions: Dict[int, str]  # position -> inserted_sequence (in donor, not in ref)
    deletions: Dict[int, str]  # position -> deleted_sequence (in ref, not in donor)


@dataclass
class EditDistanceClassification:
    """Classification result using edit distance approach."""
    outcome: str  # 'HDR_COMPLETE', 'HDR_PARTIAL', 'WT', 'NHEJ', 'MIXED'
    edits: List[DetectedEdit]
    n_donor_encoded: int  # Number of donor-encoded edits
    n_non_donor: int  # Number of non-donor edits (errors/NHEJ)
    n_total_edits: int
    donor_fraction: float  # Fraction of edits that are donor-encoded

    # For partial HDR analysis
    donor_snvs_detected: List[int] = field(default_factory=list)  # Positions of donor SNVs found
    donor_snvs_missing: List[int] = field(default_factory=list)  # Positions of donor SNVs not found


def align_donor_to_reference(
    ref_seq: str,
    donor_seq: str,
    allow_rc: bool = True
) -> DonorSignature:
    """
    Align donor to reference and extract complete edit signature (SNVs + indels).

    Uses Smith-Waterman LOCAL alignment to handle donor templates with
    5'/3' non-homologous overhangs (e.g., ssODN targeting motifs).
    This automatically finds the best matching homology region and ignores
    terminal mismatches.

    Returns:
        DonorSignature with SNVs, insertions, and deletions
    """
    ref_upper = ref_seq.upper()
    donor_upper = donor_seq.upper()

    def rev_comp(s):
        comp = {'A':'T','T':'A','G':'C','C':'G'}
        return ''.join(comp.get(b, b) for b in reversed(s.upper()))

    # Try both forward and reverse complement, use best alignment
    aligned_ref_fwd, aligned_donor_fwd, score_fwd, ref_start_fwd, donor_start_fwd = smith_waterman(ref_upper, donor_upper)

    if allow_rc:
        donor_rc = rev_comp(donor_seq)
        aligned_ref_rc, aligned_donor_rc, score_rc, ref_start_rc, donor_start_rc = smith_waterman(ref_upper, donor_rc.upper())
        if score_rc > score_fwd:
            aligned_ref, aligned_donor = aligned_ref_rc, aligned_donor_rc
            ref_start = ref_start_rc
        else:
            aligned_ref, aligned_donor = aligned_ref_fwd, aligned_donor_fwd
            ref_start = ref_start_fwd
    else:
        aligned_ref, aligned_donor = aligned_ref_fwd, aligned_donor_fwd
        ref_start = ref_start_fwd

    # Extract edits from alignment
    snvs = {}
    insertions = {}
    deletions = {}

    ref_pos = ref_start  # Position in original reference (ungapped) - START FROM ALIGNMENT OFFSET
    i = 0

    while i < len(aligned_ref):
        ref_base = aligned_ref[i]
        donor_base = aligned_donor[i]

        if ref_base == '-':
            # Insertion in donor (not in reference)
            # Collect consecutive inserted bases
            ins_seq = donor_base
            i += 1
            while i < len(aligned_ref) and aligned_ref[i] == '-':
                ins_seq += aligned_donor[i]
                i += 1
            insertions[ref_pos] = ins_seq
            continue

        if donor_base == '-':
            # Deletion in donor (present in reference)
            # Collect consecutive deleted bases
            del_start_pos = ref_pos
            del_seq = ref_base
            i += 1
            ref_pos += 1
            while i < len(aligned_ref) and aligned_donor[i] == '-':
                del_seq += aligned_ref[i]
                i += 1
                ref_pos += 1
            deletions[del_start_pos] = del_seq
            continue

        # Both present - check for SNV
        if ref_base != donor_base:
            snvs[ref_pos] = (ref_base, donor_base)

        ref_pos += 1
        i += 1

    logger.info(f"Donor alignment: {len(snvs)} SNVs, {len(insertions)} insertions, {len(deletions)} deletions at ref positions starting from {ref_start}")

    return DonorSignature(snvs=snvs, insertions=insertions, deletions=deletions)


def consolidate_donor_insertions(
    signature: 'DonorSignature',
    ref_seq: str,
    max_gap: int = 10
) -> 'DonorSignature':
    """
    Consolidate fragmented insertions into a single insertion.

    Smith-Waterman alignment sometimes produces fragmented insertions when there
    are incidental matches within a biological insertion. For example, a 10bp
    insertion like 'CGTTTCAGCT' might be split into multiple smaller insertions
    if some bases match the reference.

    This function merges insertions that are within max_gap bases of each other,
    capturing any reference bases between them as part of the consolidated insertion.

    Args:
        signature: DonorSignature from align_donor_to_reference
        ref_seq: Reference sequence (needed to fill gaps)
        max_gap: Maximum gap between insertions to merge (default 10)

    Returns:
        DonorSignature with consolidated insertions
    """
    if len(signature.insertions) <= 1:
        return signature

    # Sort insertions by position
    sorted_positions = sorted(signature.insertions.keys())

    # Find clusters of insertions to merge
    clusters = []
    current_cluster = [sorted_positions[0]]

    for pos in sorted_positions[1:]:
        prev_pos = current_cluster[-1]
        # Calculate gap (accounting for the inserted sequence length doesn't matter
        # since insertions don't consume reference positions)
        gap = pos - prev_pos

        if gap <= max_gap:
            current_cluster.append(pos)
        else:
            clusters.append(current_cluster)
            current_cluster = [pos]
    clusters.append(current_cluster)

    # Build new insertions dict
    new_insertions = {}
    ref_upper = ref_seq.upper()

    for cluster in clusters:
        if len(cluster) == 1:
            # No merging needed
            pos = cluster[0]
            new_insertions[pos] = signature.insertions[pos]
        else:
            # Merge insertions in this cluster
            # The merged insertion starts at the first position
            start_pos = cluster[0]
            _end_pos = cluster[-1]  # noqa: F841 - kept for clarity

            # Build merged sequence: insertion + ref_bases + insertion + ref_bases + ...
            merged_seq = ""
            for i, pos in enumerate(cluster):
                merged_seq += signature.insertions[pos]
                # Add reference bases between this and next insertion
                if i < len(cluster) - 1:
                    next_pos = cluster[i + 1]
                    # Reference bases between pos and next_pos
                    if pos < len(ref_upper) and next_pos <= len(ref_upper):
                        merged_seq += ref_upper[pos:next_pos]

            new_insertions[start_pos] = merged_seq
            logger.info(f"Consolidated {len(cluster)} insertions at positions {cluster} into single {len(merged_seq)}bp insertion at position {start_pos}")

    return DonorSignature(
        snvs=signature.snvs,
        insertions=new_insertions,
        deletions=signature.deletions
    )


def build_donor_signature(
    ref_seq: str,
    donor_seq: str,
    cut_site: int,
    edit_region: Optional[Tuple[int, int]] = None,
    consolidate_insertions: bool = True
) -> DonorSignature:
    """
    Build a signature of donor-encoded edits (SNVs and indels).

    Args:
        ref_seq: Reference sequence
        donor_seq: Donor template sequence
        cut_site: Position of cut site
        edit_region: Optional (start, end) to limit to specific region
        consolidate_insertions: If True, merge fragmented insertions (default True)

    Returns:
        DonorSignature containing SNVs, insertions, and deletions
    """
    signature = align_donor_to_reference(ref_seq, donor_seq)

    # Consolidate fragmented insertions (important for insertion-based edits)
    if consolidate_insertions and len(signature.insertions) > 1:
        signature = consolidate_donor_insertions(signature, ref_seq)

    # Filter to edit region if specified
    if edit_region:
        filtered_snvs = {k: v for k, v in signature.snvs.items()
                        if edit_region[0] <= k < edit_region[1]}
        filtered_ins = {k: v for k, v in signature.insertions.items()
                       if edit_region[0] <= k < edit_region[1]}
        filtered_del = {k: v for k, v in signature.deletions.items()
                       if edit_region[0] <= k < edit_region[1]}
        return DonorSignature(snvs=filtered_snvs, insertions=filtered_ins, deletions=filtered_del)

    return signature


def _detect_homopolymer_length(ref_seq: str, pos: int) -> int:
    """Detect length of homopolymer run at given position."""
    if pos < 0 or pos >= len(ref_seq):
        return 0
    base = ref_seq[pos].upper()
    run_len = 1
    # Extend forward
    for i in range(pos + 1, min(pos + 15, len(ref_seq))):
        if ref_seq[i].upper() == base:
            run_len += 1
        else:
            break
    # Extend backward
    for i in range(pos - 1, max(pos - 15, -1), -1):
        if ref_seq[i].upper() == base:
            run_len += 1
        else:
            break
    return run_len


def _get_indel_microhomology_lengths(
    indel_pos: int,
    indel_length: int,
    indel_type: EditType,
    indel_seq: str,
    ref_seq: str
) -> Tuple[int, int]:
    """
    Get microhomology lengths at indel boundaries.

    Returns the length of microhomology on the left and right sides of an indel.
    Used to distinguish MMEJ (>2bp microhomology) from classical NHEJ (0-2bp).

    Args:
        indel_pos: Reported position of indel in reference
        indel_length: Length of indel
        indel_type: EditType.INSERTION or EditType.DELETION
        indel_seq: Sequence of indel (deleted or inserted bases)
        ref_seq: Reference sequence

    Returns:
        (left_microhomology_length, right_microhomology_length) tuple
    """
    ref_upper = ref_seq.upper()
    indel_seq_upper = indel_seq.upper()

    if indel_type == EditType.DELETION:
        deleted_seq = ref_upper[indel_pos:indel_pos + indel_length] if indel_pos + indel_length <= len(ref_upper) else ""
        if not deleted_seq:
            return (0, 0)

        # Left microhomology: how far left could deletion be shifted
        left_shift = 0
        for i in range(1, min(indel_length, indel_pos) + 1):
            if indel_pos - i >= 0 and deleted_seq[:i] == ref_upper[indel_pos - i:indel_pos]:
                left_shift = i

        # Right microhomology: how far right could deletion be shifted
        right_shift = 0
        for i in range(1, min(indel_length, len(ref_upper) - indel_pos - indel_length) + 1):
            if indel_pos + indel_length + i <= len(ref_upper) and deleted_seq[-i:] == ref_upper[indel_pos + indel_length:indel_pos + indel_length + i]:
                right_shift = i

        return (left_shift, right_shift)

    else:  # INSERTION
        if not indel_seq_upper:
            return (0, 0)

        ins_len = len(indel_seq_upper)

        # Left microhomology: how far left could insertion be shifted
        left_extend = 0
        for i in range(1, min(ins_len, indel_pos) + 1):
            if indel_pos - i >= 0 and indel_seq_upper[-i:] == ref_upper[indel_pos - i:indel_pos]:
                left_extend = i

        # Right microhomology: how far right could insertion be shifted
        right_extend = 0
        for i in range(1, min(ins_len, len(ref_upper) - indel_pos) + 1):
            if indel_pos + i <= len(ref_upper) and indel_seq_upper[:i] == ref_upper[indel_pos:indel_pos + i]:
                right_extend = i

        return (left_extend, right_extend)


def _get_indel_position_range(
    indel_pos: int,
    indel_length: int,
    indel_type: EditType,
    indel_seq: str,
    ref_seq: str
) -> Tuple[int, int]:
    """
    Get the range of positions where this indel could be placed due to microhomology.

    Accounts for alignment ambiguity in repetitive regions. The same biological
    deletion can be left-aligned, right-aligned, or anywhere in between if there
    is microhomology at the boundaries.

    Args:
        indel_pos: Reported position of indel in reference
        indel_length: Length of indel
        indel_type: EditType.INSERTION or EditType.DELETION
        indel_seq: Sequence of indel (deleted or inserted bases)
        ref_seq: Reference sequence

    Returns:
        (min_pos, max_pos) tuple representing the range where this indel could occur
    """
    left_shift, right_shift = _get_indel_microhomology_lengths(
        indel_pos, indel_length, indel_type, indel_seq, ref_seq
    )
    ref_upper = ref_seq.upper()
    indel_seq_upper = indel_seq.upper()

    if indel_type == EditType.DELETION:
        # For deletions: check microhomology at boundaries
        deleted_seq = ref_upper[indel_pos:indel_pos + indel_length] if indel_pos + indel_length <= len(ref_upper) else ""

        if not deleted_seq:
            return (indel_pos, indel_pos)

        # How far left could this deletion be shifted?
        # Check ALL positions (don't break on first mismatch) to handle interrupted repeats
        left_shift = 0
        for i in range(1, min(indel_length, indel_pos) + 1):
            # Check if prefix of deleted seq matches suffix before deletion
            if indel_pos - i >= 0 and deleted_seq[:i] == ref_upper[indel_pos - i:indel_pos]:
                left_shift = i
            # NO BREAK - continue checking to find maximum shift (handles interrupted homopolymers)

        # How far right could this deletion be shifted?
        right_shift = 0
        for i in range(1, min(indel_length, len(ref_upper) - indel_pos - indel_length) + 1):
            # Check if suffix of deleted seq matches prefix after deletion
            if indel_pos + indel_length + i <= len(ref_upper) and deleted_seq[-i:] == ref_upper[indel_pos + indel_length:indel_pos + indel_length + i]:
                right_shift = i
            # NO BREAK - continue checking to find maximum shift

        # Deletion spans [indel_pos, indel_pos + length)
        # It could be shifted left by left_shift or right by right_shift
        return (indel_pos - left_shift, indel_pos + indel_length + right_shift)

    else:  # INSERTION
        # Check if insertion sequence matches tandem repeats in the reference
        # This handles homopolymers (AAAAA), dinucleotide repeats (ATATAT),
        # trinucleotide repeats (ATCATCATC), etc.
        if not indel_seq_upper:
            return (indel_pos, indel_pos)

        ins_len = len(indel_seq_upper)

        # Check how far left the insertion could be shifted
        # It can shift if the reference sequence before insertion matches the inserted sequence
        left_extend = 0
        for shift in range(1, min(ins_len * 5, indel_pos) + 1):  # Check up to 5x insertion length
            # Check if we can shift by 'shift' bases
            # This requires that ref[pos-shift:pos] contains the same pattern as the insertion
            check_start = indel_pos - shift
            check_end = indel_pos
            if check_start >= 0:
                ref_context = ref_upper[check_start:check_end]
                # Check if this context is compatible with the insertion (tandem repeat)
                if ref_context == indel_seq_upper[-shift:]:
                    left_extend = shift
                # NO BREAK - continue to find maximum shift (handles interrupted repeats)

        # Check how far right the insertion could be shifted
        right_extend = 0
        for shift in range(1, min(ins_len * 5, len(ref_upper) - indel_pos) + 1):
            check_start = indel_pos
            check_end = indel_pos + shift
            if check_end <= len(ref_upper):
                ref_context = ref_upper[check_start:check_end]
                # Check if this context is compatible with the insertion (tandem repeat)
                if ref_context == indel_seq_upper[:shift]:
                    right_extend = shift
                # NO BREAK - continue to find maximum shift

        # Insertion could slide within the tandem repeat region
        return (indel_pos - left_extend, indel_pos + right_extend)


def _could_involve_cut_site(
    indel_pos: int,
    indel_length: int,
    indel_type: EditType,
    indel_seq: str,
    ref_seq: str,
    cut_site: int,
    window: int = 1
) -> bool:
    """
    Check if indel could plausibly result from NHEJ repair at cut_site.

    NHEJ indels MUST involve the cut site. The cut site is at a specific position
    (e.g., 222). An indel can only be classified as NHEJ if its microhomology range
    includes that position.

    The microhomology analysis accounts for tandem repeats:
    - A 1bp deletion at pos 268 in a 7T run → range [262, 274]
    - Cut site at 222 is NOT in [262, 274] → NOT NHEJ
    - A 1bp deletion at pos 220 in a 2T run → range [220, 221]
    - Cut site at 222 is NOT in [220, 221] → NOT NHEJ (too far!)
    - A 2bp deletion at pos 221 → range [221, 223]
    - Cut site at 222 IS in [221, 223] → IS NHEJ!

    Args:
        indel_pos: Reported position of indel
        indel_length: Length of indel
        indel_type: EditType.INSERTION or EditType.DELETION
        indel_seq: Sequence of indel
        ref_seq: Reference sequence
        cut_site: Position of Cas9 cut site (single coordinate)
        window: Tiny window for Cas9 cutting position uncertainty (default: 1bp).
                SpCas9 cuts 3bp upstream of PAM, but position can vary ±1bp.

    Returns:
        True if indel's microhomology range includes the cut site
    """
    # Get the range where this indel could be placed due to tandem repeats
    # This ALREADY accounts for microhomology-mediated positional ambiguity
    min_pos, max_pos = _get_indel_position_range(
        indel_pos, indel_length, indel_type, indel_seq, ref_seq
    )

    # Check if the cut site falls within the indel's microhomology range
    # Allow a tiny window (±1bp) for Cas9 cutting position variability
    return (min_pos <= cut_site + window) and (max_pos >= cut_site - window)


def classify_read_edit_distance(
    read_seq: str,
    ref_seq: str,
    donor_signature: DonorSignature,
    ref_start: int,
    cigar_ops: Optional[List[Tuple[int, int]]] = None,
    cut_site: int = 0,
    snv_distance_filter: int = 50,
    filter_chimeric: bool = True,
    homopolymer_filter: int = 0,
    nhej_quantification_window: int = 1,
    min_alignment_quality: float = 0.90,
    max_mismatch_rate: float = 0.15
) -> EditDistanceClassification:
    """
    Classify a read using edit distance approach.

    For each difference between read and reference:
    - Check if it matches the donor signature (SNVs or indels)
    - If yes, it's donor-encoded
    - If no, classify based on edit type and location

    Outcome categories:
    - WT: no edits detected (includes perfect NHEJ repairs)
    - HDR_COMPLETE: all donor SNVs/indels present, no other edits
    - HDR_PARTIAL: some donor SNVs/indels present
    - NHEJ_INDEL: classical NHEJ indels (0-2bp microhomology) at cut site, no donor
    - MMEJ_INDEL: alt-NHEJ indels (>2bp microhomology) at cut site, no donor
    - HDR_PLUS_NHEJ_INDEL: donor markers + classical NHEJ indels at cut site
    - HDR_PLUS_MMEJ_INDEL: donor markers + MMEJ indels at cut site
    - HDR_PLUS_OTHER: donor markers + edits NOT at cut site
    - NON_DONOR_SNV: SNVs not in donor, no indels
    - NON_DONOR_NON_NHEJ_INDEL: indels NOT at cut site
    - LOW_QUALITY: alignment quality below threshold (<90% match or >15% mismatch),
                   likely sequencing errors or poor read quality
    - CHIMERIC: low-quality artifacts (large insertions with N runs, likely
                sequencing quality issues, adapter contamination, or PCR artifacts)

    Microhomology-based pathway classification:
    - Classical NHEJ: 0-2bp microhomology (canonical DNA-PKcs/Ligase IV pathway)
    - MMEJ (Alt-NHEJ): >2bp microhomology (Polθ-mediated alternative pathway)

    Args:
        read_seq: Read sequence
        ref_seq: Reference sequence
        donor_signature: DonorSignature with SNVs, insertions, and deletions
        ref_start: Start position of read in reference
        cigar_ops: Optional CIGAR operations for indel handling
        cut_site: Cut site position for distance calculations
        snv_distance_filter: Only count non-donor SNVs within this distance
            of cut_site as potentially meaningful (default: 50bp). SNVs farther
            away are likely sequencing errors and ignored for classification.
        filter_chimeric: If True, reads with large insertions (>=30bp) containing
            N runs (>=5 N's) are classified as CHIMERIC artifacts and excluded from
            analysis. N's typically indicate low base quality or sequencing artifacts
            rather than genuine biological variation (default: True).
        homopolymer_filter: If >0, 1bp indels at homopolymer runs of this length
            or greater are ignored as sequencing artifacts (default: 0 = disabled).
            Recommended: 5-7 for Illumina data.
        nhej_quantification_window: Window for Cas9 cutting position uncertainty
            (default: 1bp). SpCas9 cuts 3bp upstream of PAM, but position can vary
            ±1bp. Indels must have microhomology range that includes cut_site ± window.
        min_alignment_quality: Minimum fraction of matching bases required for valid
            alignment (default: 0.90 = 90% match). Reads below this threshold are
            classified as LOW_QUALITY to prevent false positive HDR calls on poor
            alignments. Quality = (total_length - mismatches) / total_length.
        max_mismatch_rate: Maximum fraction of mismatched bases allowed (default: 0.15
            = 15% mismatch). Alternative threshold to min_alignment_quality. Reads
            exceeding this are classified as LOW_QUALITY.

    Returns:
        EditDistanceClassification
    """
    read_upper = read_seq.upper()
    ref_upper = ref_seq.upper()

    edits = []
    donor_snvs_detected = []
    donor_snvs_missing = list(donor_signature.snvs.keys())
    is_chimeric = False  # Flag for chimeric read detection

    # Process CIGAR if available (for proper indel handling)
    if cigar_ops:
        ref_pos = ref_start
        read_pos = 0

        for op, length in cigar_ops:
            if op == 0:  # M (match/mismatch)
                for i in range(length):
                    if read_pos < len(read_upper) and ref_pos < len(ref_upper):
                        read_base = read_upper[read_pos]
                        ref_base = ref_upper[ref_pos]

                        if read_base != ref_base:
                            # Check if this is donor-encoded SNV
                            is_donor = False
                            if ref_pos in donor_signature.snvs:
                                expected_ref, expected_donor = donor_signature.snvs[ref_pos]
                                if read_base == expected_donor:
                                    is_donor = True
                                    donor_snvs_detected.append(ref_pos)
                                    if ref_pos in donor_snvs_missing:
                                        donor_snvs_missing.remove(ref_pos)

                            edits.append(DetectedEdit(
                                edit_type=EditType.SNV,
                                ref_position=ref_pos,
                                ref_base=ref_base,
                                read_base=read_base,
                                is_donor_encoded=is_donor,
                                distance_to_cut=ref_pos - cut_site
                            ))

                    ref_pos += 1
                    read_pos += 1

            elif op == 1:  # I (insertion)
                inserted_seq = read_upper[read_pos:read_pos + length]

                # Check for low-quality artifact: large insertion with N runs
                # These are likely: (1) low base quality calls, (2) adapter contamination,
                # (3) PCR artifacts, or (4) index hopping. N's typically indicate
                # sequencing quality issues rather than genuine biological sequences.
                if filter_chimeric and length >= 30 and inserted_seq.count('N') >= 5:
                    is_chimeric = True

                # Skip 1bp insertions at homopolymer positions (sequencing artifact)
                skip_homopolymer = False
                if homopolymer_filter > 0 and length == 1:
                    hp_len = _detect_homopolymer_length(ref_upper, ref_pos)
                    if hp_len >= homopolymer_filter:
                        skip_homopolymer = True

                if not skip_homopolymer:
                    # Check if insertion matches donor
                    # Use flexible matching: position within window AND sequence overlap
                    is_donor_ins = False
                    position_tolerance = 5  # Allow insertions within ±5bp of expected

                    for expected_pos, expected_ins in donor_signature.insertions.items():
                        # Check if position is close
                        if abs(ref_pos - expected_pos) <= position_tolerance:
                            # Check for sequence match (exact, substring, or high overlap)
                            if inserted_seq == expected_ins:
                                # Exact match
                                is_donor_ins = True
                                break
                            elif inserted_seq in expected_ins:
                                # Read insertion is subset of expected (partial HDR)
                                is_donor_ins = True
                                break
                            elif expected_ins in inserted_seq:
                                # Read insertion contains expected (HDR with extra bases)
                                is_donor_ins = True
                                break
                            else:
                                # Check overlap for long insertions (>5bp)
                                if len(inserted_seq) >= 5 and len(expected_ins) >= 5:
                                    # Check if core sequence matches
                                    min_overlap = min(len(inserted_seq), len(expected_ins)) // 2
                                    for offset in range(-3, 4):
                                        if offset >= 0:
                                            ins_slice = inserted_seq[offset:]
                                            exp_slice = expected_ins[:len(ins_slice)]
                                        else:
                                            ins_slice = inserted_seq[:len(expected_ins)+offset]
                                            exp_slice = expected_ins[-offset:]
                                        if len(ins_slice) >= min_overlap and ins_slice == exp_slice:
                                            is_donor_ins = True
                                            break
                                    if is_donor_ins:
                                        break

                    edits.append(DetectedEdit(
                        edit_type=EditType.INSERTION,
                        ref_position=ref_pos,
                        ref_base="",
                        read_base=inserted_seq,
                        is_donor_encoded=is_donor_ins,
                        distance_to_cut=ref_pos - cut_site
                    ))
                read_pos += length

            elif op == 2:  # D (deletion)
                deleted_seq = ref_upper[ref_pos:ref_pos + length]

                # Skip 1bp deletions at homopolymer positions (sequencing artifact)
                skip_homopolymer = False
                if homopolymer_filter > 0 and length == 1:
                    hp_len = _detect_homopolymer_length(ref_upper, ref_pos)
                    if hp_len >= homopolymer_filter:
                        skip_homopolymer = True

                if not skip_homopolymer:
                    # Check if deletion matches donor
                    is_donor_del = False
                    if ref_pos in donor_signature.deletions:
                        expected_del = donor_signature.deletions[ref_pos]
                        if deleted_seq == expected_del:
                            is_donor_del = True

                    edits.append(DetectedEdit(
                        edit_type=EditType.DELETION,
                        ref_position=ref_pos,
                        ref_base=deleted_seq,
                        read_base="",
                        is_donor_encoded=is_donor_del,
                        distance_to_cut=ref_pos - cut_site
                    ))
                ref_pos += length

            elif op == 4:  # S (soft clip)
                read_pos += length
            elif op == 5:  # H (hard clip)
                pass
    else:
        # Simple comparison without CIGAR (assumes no indels)
        for i, read_base in enumerate(read_upper):
            ref_pos = ref_start + i
            if ref_pos >= len(ref_upper):
                break

            ref_base = ref_upper[ref_pos]

            if read_base != ref_base:
                # Check if donor-encoded
                is_donor = False
                if ref_pos in donor_signature:
                    expected_ref, expected_donor = donor_signature[ref_pos]
                    if read_base == expected_donor:
                        is_donor = True
                        donor_snvs_detected.append(ref_pos)
                        if ref_pos in donor_snvs_missing:
                            donor_snvs_missing.remove(ref_pos)

                edits.append(DetectedEdit(
                    edit_type=EditType.SNV,
                    ref_position=ref_pos,
                    ref_base=ref_base,
                    read_base=read_base,
                    is_donor_encoded=is_donor,
                    distance_to_cut=ref_pos - cut_site
                ))

    # Calculate alignment quality to filter low-quality reads
    # Total mismatches = SNVs + insertion bases + deletion bases
    total_snvs = sum(1 for e in edits if e.edit_type == EditType.SNV)
    total_insertion_bases = sum(len(e.read_base) for e in edits if e.edit_type == EditType.INSERTION)
    total_deletion_bases = sum(len(e.ref_base) for e in edits if e.edit_type == EditType.DELETION)
    total_mismatches = total_snvs + total_insertion_bases + total_deletion_bases

    # Alignment length is the longer of read or reference sequence (accounting for indels)
    alignment_length = max(len(read_upper), len(ref_upper))

    # Calculate quality metrics
    if alignment_length > 0:
        alignment_quality = (alignment_length - total_mismatches) / alignment_length
        mismatch_rate = total_mismatches / alignment_length
    else:
        alignment_quality = 0.0
        mismatch_rate = 1.0

    # Filter out low-quality alignments that could produce false positive HDR calls
    if alignment_quality < min_alignment_quality or mismatch_rate > max_mismatch_rate:
        logger.debug(f"LOW_QUALITY alignment: quality={alignment_quality:.3f}, mismatch_rate={mismatch_rate:.3f}, "
                    f"mismatches={total_mismatches}, length={alignment_length}")
        return EditDistanceClassification(
            outcome='LOW_QUALITY',
            edits=edits,
            n_donor_encoded=0,
            n_non_donor=len(edits),
            n_total_edits=len(edits),
            donor_fraction=0.0,
            donor_snvs_detected=[],
            donor_snvs_missing=list(donor_signature.snvs.keys())
        )

    # Check for chimeric reads first (index hopping artifacts)
    if is_chimeric:
        return EditDistanceClassification(
            outcome='CHIMERIC',
            edits=edits,
            n_donor_encoded=0,
            n_non_donor=len(edits),
            n_total_edits=len(edits),
            donor_fraction=0.0,
            donor_snvs_detected=[],
            donor_snvs_missing=list(donor_signature.snvs.keys())
        )

    # Count donor vs non-donor edits
    n_donor = sum(1 for e in edits if e.is_donor_encoded)
    n_non_donor = len(edits) - n_donor
    donor_fraction = n_donor / len(edits) if edits else 0

    # Separate indels by whether they could involve the cut site (NHEJ) or not
    nhej_indels = []
    non_nhej_indels = []

    for e in edits:
        if e.edit_type != EditType.SNV and not e.is_donor_encoded:
            # Check if this indel could plausibly involve the cut site
            indel_seq = e.read_base if e.edit_type == EditType.INSERTION else e.ref_base
            could_be_nhej = _could_involve_cut_site(
                indel_pos=e.ref_position,
                indel_length=e.size,
                indel_type=e.edit_type,
                indel_seq=indel_seq,
                ref_seq=ref_upper,
                cut_site=cut_site,
                window=nhej_quantification_window
            )
            if could_be_nhej:
                nhej_indels.append(e)
            else:
                non_nhej_indels.append(e)

    # Further classify NHEJ indels into MMEJ vs classical NHEJ based on microhomology
    # MMEJ: microhomology > 2bp on EITHER side (indicates alt-NHEJ pathway)
    # Classical NHEJ: microhomology 0-2bp (indicates canonical NHEJ pathway)
    mmej_indels = []
    classical_nhej_indels = []

    for e in nhej_indels:
        indel_seq = e.read_base if e.edit_type == EditType.INSERTION else e.ref_base
        left_mh, right_mh = _get_indel_microhomology_lengths(
            e.ref_position, e.size, e.edit_type, indel_seq, ref_upper
        )
        # MMEJ: microhomology > 2bp on EITHER side
        if max(left_mh, right_mh) > 2:
            mmej_indels.append(e)
        else:
            classical_nhej_indels.append(e)

    n_mmej_indels = len(mmej_indels)
    n_classical_nhej_indels = len(classical_nhej_indels)
    n_non_nhej_indels = len(non_nhej_indels)

    # Only count non-donor SNVs within snv_distance_filter of cut site
    # SNVs farther away are likely sequencing errors, not editing outcomes
    n_non_donor_snvs = sum(
        1 for e in edits
        if e.edit_type == EditType.SNV
        and not e.is_donor_encoded
        and abs(e.distance_to_cut) <= snv_distance_filter
    )

    # Determine outcome
    # Count all expected donor edits (SNVs + indels)
    n_expected_donor_edits = (len(donor_signature.snvs) +
                              len(donor_signature.insertions) +
                              len(donor_signature.deletions))

    if len(edits) == 0:
        outcome = 'WT'
    elif n_donor == n_expected_donor_edits and n_non_donor == 0:
        outcome = 'HDR_COMPLETE'
    elif n_donor > 0 and n_non_donor == 0:
        # Has some donor edits (SNVs or indels), no other edits
        outcome = 'HDR_PARTIAL'
    elif n_donor > 0 and n_mmej_indels > 0:
        # Donor SNVs + MMEJ indels at cut site (prioritize MMEJ)
        outcome = 'HDR_PLUS_MMEJ_INDEL'
    elif n_donor > 0 and n_classical_nhej_indels > 0:
        # Donor SNVs + classical NHEJ indels at cut site
        outcome = 'HDR_PLUS_NHEJ_INDEL'
    elif n_donor > 0 and (n_non_nhej_indels > 0 or n_non_donor_snvs > 0):
        # Donor SNVs + edits NOT at cut site
        outcome = 'HDR_PLUS_OTHER'
    elif n_mmej_indels > 0 and n_donor == 0:
        # MMEJ indels at cut site, no donor markers (alt-NHEJ repair pathway)
        outcome = 'MMEJ_INDEL'
    elif n_classical_nhej_indels > 0 and n_donor == 0:
        # Classical NHEJ indels at cut site, no donor markers (canonical NHEJ repair)
        outcome = 'NHEJ_INDEL'
    elif n_non_nhej_indels > 0 and n_donor == 0:
        # Indels NOT at cut site, no donor markers
        outcome = 'NON_DONOR_NON_NHEJ_INDEL'
    elif n_donor == 0 and n_non_donor_snvs > 0:
        # Only has SNVs (no donor, no indels)
        outcome = 'NON_DONOR_SNV'
    else:
        # No edits or only distant SNVs (treated as sequencing errors)
        outcome = 'WT'

    return EditDistanceClassification(
        outcome=outcome,
        edits=edits,
        n_donor_encoded=n_donor,
        n_non_donor=n_non_donor,
        n_total_edits=len(edits),
        donor_fraction=donor_fraction,
        donor_snvs_detected=donor_snvs_detected,
        donor_snvs_missing=donor_snvs_missing
    )


def analyze_donor_integration(
    classifications: List[EditDistanceClassification],
    donor_signature: Dict[int, Tuple[str, str]],
    cut_site: int
) -> Dict:
    """
    Analyze donor integration patterns across all reads.

    Returns:
        Dict with:
        - per_position_frequency: Dict[position -> frequency]
        - distance_gradient: frequency vs distance from cut
        - partial_hdr_stats: breakdown of partial HDR
    """
    n_reads = len(classifications)

    # Count per-position integration
    position_counts = {pos: 0 for pos in donor_signature.keys()}

    hdr_reads = 0
    for clf in classifications:
        if clf.outcome in ('HDR_COMPLETE', 'HDR_PARTIAL', 'MIXED'):
            hdr_reads += 1
            for pos in clf.donor_snvs_detected:
                if pos in position_counts:
                    position_counts[pos] += 1

    # Calculate frequencies
    per_position_freq = {}
    for pos, (ref_base, donor_base) in donor_signature.items():
        count = position_counts[pos]
        freq = count / hdr_reads if hdr_reads > 0 else 0
        distance = pos - cut_site
        per_position_freq[pos] = {
            'frequency': freq,
            'count': count,
            'distance_to_cut': distance,
            'ref_base': ref_base,
            'donor_base': donor_base
        }

    # Partial HDR breakdown
    partial_counts = {}
    for clf in classifications:
        if clf.outcome == 'HDR_PARTIAL':
            n = len(clf.donor_snvs_detected)
            partial_counts[n] = partial_counts.get(n, 0) + 1

    # Outcome summary
    outcomes = {}
    for clf in classifications:
        outcomes[clf.outcome] = outcomes.get(clf.outcome, 0) + 1

    return {
        'total_reads': n_reads,
        'hdr_reads': hdr_reads,
        'outcomes': outcomes,
        'per_position_frequency': per_position_freq,
        'partial_breakdown': partial_counts,
        'n_donor_positions': len(donor_signature)
    }


@dataclass
class MultiDonorClassification:
    """
    Classification result when comparing against multiple donor templates.

    Used for barcoded experiments where each sample has a unique HDR template
    and we need to:
    1. Find the best matching donor from all possibilities
    2. Detect sample swaps (when detected donor differs from expected)
    """
    # Best matching donor
    best_donor_id: str
    best_donor_score: float  # Higher is better match
    best_donor_classification: EditDistanceClassification

    # Expected donor (from manifest)
    expected_donor_id: str

    # Sample swap detection
    is_sample_swap: bool  # True if best_donor_id != expected_donor_id
    swap_confidence: float  # Confidence in swap detection (0-1)

    # All donor scores for debugging/analysis
    all_donor_scores: Dict[str, float] = field(default_factory=dict)


def score_donor_match(
    classification: EditDistanceClassification,
    donor_signature: DonorSignature
) -> float:
    """
    Score how well a read matches a donor template.

    Scoring system:
    - +10 points per donor-encoded SNV detected
    - +10 points per donor-encoded insertion detected
    - +10 points per donor-encoded deletion detected
    - -5 points per non-donor edit
    - Bonus for complete HDR: +20 points

    Args:
        classification: Result from classify_read_edit_distance
        donor_signature: The donor template signature

    Returns:
        Score (higher = better match)
    """
    score = 0.0

    # Count expected donor edits
    n_expected_snvs = len(donor_signature.snvs)
    n_expected_insertions = len(donor_signature.insertions)
    n_expected_deletions = len(donor_signature.deletions)
    n_total_expected = n_expected_snvs + n_expected_insertions + n_expected_deletions

    if n_total_expected == 0:
        # No donor edits expected (WT/control)
        return 0.0

    # Count donor edits found
    n_donor_edits = classification.n_donor_encoded
    n_non_donor = classification.n_non_donor

    # Base score: points for donor edits found
    score += n_donor_edits * 10.0

    # Penalty for non-donor edits
    score -= n_non_donor * 5.0

    # Bonus for complete HDR
    if n_donor_edits == n_total_expected and n_non_donor == 0:
        score += 20.0

    # Normalize by expected edits to handle different barcode sizes
    # This makes scores comparable across donors with different numbers of edits
    normalized_score = score / n_total_expected

    return normalized_score


def classify_read_multi_donor(
    read_seq: str,
    ref_seq: str,
    donor_signatures: Dict[str, DonorSignature],
    expected_donor_id: str,
    ref_start: int,
    cigar_ops: Optional[List[Tuple[int, int]]] = None,
    cut_site: int = 0,
    snv_distance_filter: int = 50,
    filter_chimeric: bool = True,
    homopolymer_filter: int = 0,
    nhej_quantification_window: int = 1,
    min_alignment_quality: float = 0.90,
    max_mismatch_rate: float = 0.15,
    min_score_difference: float = 5.0
) -> MultiDonorClassification:
    """
    Classify a read by comparing against multiple donor templates.

    This enables:
    1. Sample swap detection: Find which barcode is actually present
    2. Improved HDR calling: Match against the correct template
    3. Quality control: Flag ambiguous assignments

    Args:
        read_seq: Read sequence
        ref_seq: Reference sequence
        donor_signatures: Dict mapping donor_id -> DonorSignature
        expected_donor_id: The donor ID expected for this sample (from manifest)
        ref_start: Start position of read in reference
        cigar_ops: Optional CIGAR operations for indel handling
        cut_site: Cut site position for distance calculations
        snv_distance_filter: Distance filter for non-donor SNVs
        filter_chimeric: Filter chimeric reads
        homopolymer_filter: Filter homopolymer indels
        nhej_quantification_window: Window for NHEJ calling
        min_alignment_quality: Minimum alignment quality
        max_mismatch_rate: Maximum mismatch rate
        min_score_difference: Minimum score difference to confidently call a swap

    Returns:
        MultiDonorClassification with best donor match and swap detection
    """
    if not donor_signatures:
        # No donors to compare - return empty result
        empty_clf = EditDistanceClassification(
            outcome='NO_DONORS',
            edits=[],
            n_donor_encoded=0,
            n_non_donor=0,
            n_total_edits=0,
            donor_fraction=0.0
        )
        return MultiDonorClassification(
            best_donor_id='',
            best_donor_score=0.0,
            best_donor_classification=empty_clf,
            expected_donor_id=expected_donor_id,
            is_sample_swap=False,
            swap_confidence=0.0
        )

    # Classify read against each donor
    donor_results = {}
    donor_scores = {}

    for donor_id, signature in donor_signatures.items():
        clf = classify_read_edit_distance(
            read_seq=read_seq,
            ref_seq=ref_seq,
            donor_signature=signature,
            ref_start=ref_start,
            cigar_ops=cigar_ops,
            cut_site=cut_site,
            snv_distance_filter=snv_distance_filter,
            filter_chimeric=filter_chimeric,
            homopolymer_filter=homopolymer_filter,
            nhej_quantification_window=nhej_quantification_window,
            min_alignment_quality=min_alignment_quality,
            max_mismatch_rate=max_mismatch_rate
        )

        score = score_donor_match(clf, signature)
        donor_results[donor_id] = clf
        donor_scores[donor_id] = score

    # Find best matching donor
    if donor_scores:
        best_donor_id = max(donor_scores, key=donor_scores.get)
        best_score = donor_scores[best_donor_id]
        best_clf = donor_results[best_donor_id]

        # Calculate swap confidence
        # If best donor differs from expected, calculate confidence
        if best_donor_id != expected_donor_id:
            expected_score = donor_scores.get(expected_donor_id, 0.0)
            score_diff = best_score - expected_score

            # Only call swap if score difference exceeds threshold
            is_swap = score_diff >= min_score_difference and best_score > 0

            # Confidence based on score difference
            if is_swap:
                # Normalize confidence: score_diff of min_score_difference = 0.5,
                # score_diff of 2*min_score_difference = 0.75, etc.
                swap_confidence = min(0.5 + (score_diff - min_score_difference) / (2 * min_score_difference), 1.0)
            else:
                swap_confidence = 0.0
        else:
            is_swap = False
            swap_confidence = 0.0
    else:
        # Fallback to expected donor if no scores
        best_donor_id = expected_donor_id
        best_score = 0.0
        best_clf = donor_results.get(expected_donor_id, EditDistanceClassification(
            outcome='NO_MATCH',
            edits=[],
            n_donor_encoded=0,
            n_non_donor=0,
            n_total_edits=0,
            donor_fraction=0.0
        ))
        is_swap = False
        swap_confidence = 0.0

    return MultiDonorClassification(
        best_donor_id=best_donor_id,
        best_donor_score=best_score,
        best_donor_classification=best_clf,
        expected_donor_id=expected_donor_id,
        is_sample_swap=is_swap,
        swap_confidence=swap_confidence,
        all_donor_scores=donor_scores
    )


def batch_classify_multi_donor(
    sequences: List[Tuple[str, str, int, Optional[List[Tuple[int, int]]]]],
    ref_seq: str,
    donor_signatures: Dict[str, DonorSignature],
    expected_donor_id: str,
    cut_site: int = 0,
    **kwargs
) -> List[MultiDonorClassification]:
    """
    Classify multiple sequences against multiple donors efficiently.

    Optimizes by pre-building all donor signatures once, then classifying
    each sequence against all donors.

    Args:
        sequences: List of (read_seq, seq_id, ref_start, cigar_ops) tuples
        ref_seq: Reference sequence
        donor_signatures: Dict mapping donor_id -> DonorSignature
        expected_donor_id: Expected donor for this sample
        cut_site: Cut site position
        **kwargs: Additional arguments for classify_read_multi_donor

    Returns:
        List of MultiDonorClassification results
    """
    results = []

    for read_seq, seq_id, ref_start, cigar_ops in sequences:
        clf = classify_read_multi_donor(
            read_seq=read_seq,
            ref_seq=ref_seq,
            donor_signatures=donor_signatures,
            expected_donor_id=expected_donor_id,
            ref_start=ref_start,
            cigar_ops=cigar_ops,
            cut_site=cut_site,
            **kwargs
        )
        results.append(clf)

    return results


def summarize_sample_swap_detection(
    classifications: List[MultiDonorClassification]
) -> Dict:
    """
    Summarize sample swap detection results across all reads in a sample.

    Args:
        classifications: List of MultiDonorClassification results

    Returns:
        Dict with:
        - expected_donor_id: Expected donor for this sample
        - detected_donor_id: Most commonly detected donor
        - is_sample_swap: Whether sample appears to be swapped
        - swap_confidence: Confidence in swap detection
        - donor_distribution: Dict of donor_id -> read count
        - swap_reads_fraction: Fraction of reads matching non-expected donor
    """
    if not classifications:
        return {
            'expected_donor_id': '',
            'detected_donor_id': '',
            'is_sample_swap': False,
            'swap_confidence': 0.0,
            'donor_distribution': {},
            'swap_reads_fraction': 0.0
        }

    expected_donor_id = classifications[0].expected_donor_id

    # Count reads per detected donor (only count HDR-positive reads)
    donor_counts = {}
    total_hdr_reads = 0

    for clf in classifications:
        # Only count reads that show HDR evidence
        if clf.best_donor_score > 0 and clf.best_donor_classification.n_donor_encoded > 0:
            best_id = clf.best_donor_id
            donor_counts[best_id] = donor_counts.get(best_id, 0) + 1
            total_hdr_reads += 1

    if total_hdr_reads == 0:
        return {
            'expected_donor_id': expected_donor_id,
            'detected_donor_id': expected_donor_id,
            'is_sample_swap': False,
            'swap_confidence': 0.0,
            'donor_distribution': {},
            'swap_reads_fraction': 0.0
        }

    # Find most common donor
    detected_donor_id = max(donor_counts, key=donor_counts.get)
    detected_count = donor_counts[detected_donor_id]

    # Calculate swap metrics
    expected_count = donor_counts.get(expected_donor_id, 0)
    swap_reads_fraction = 1.0 - (expected_count / total_hdr_reads)

    # Call swap if detected donor differs and has clear majority
    is_swap = (detected_donor_id != expected_donor_id and
               detected_count > expected_count and
               swap_reads_fraction > 0.5)

    # Confidence based on how dominant the detected donor is
    if is_swap:
        swap_confidence = detected_count / total_hdr_reads
    else:
        swap_confidence = 0.0

    return {
        'expected_donor_id': expected_donor_id,
        'detected_donor_id': detected_donor_id,
        'is_sample_swap': is_swap,
        'swap_confidence': swap_confidence,
        'donor_distribution': donor_counts,
        'swap_reads_fraction': swap_reads_fraction,
        'total_hdr_reads': total_hdr_reads
    }
