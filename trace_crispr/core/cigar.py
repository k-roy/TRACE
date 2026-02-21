"""
CIGAR string parsing utilities for editing outcome analysis.

Author: Kevin R. Roy
"""

from dataclasses import dataclass
from typing import List, Optional, Tuple

import pysam

# CIGAR operation codes (from pysam/SAM spec)
CIGAR_OPS = {
    0: 'M',   # Match/mismatch
    1: 'I',   # Insertion
    2: 'D',   # Deletion
    3: 'N',   # Skipped region (intron)
    4: 'S',   # Soft clip
    5: 'H',   # Hard clip
    6: 'P',   # Padding
    7: '=',   # Sequence match
    8: 'X',   # Sequence mismatch
}

# Operations that consume reference bases
REF_CONSUMING_OPS = {0, 2, 3, 7, 8}  # M, D, N, =, X

# Operations that consume query (read) bases
QUERY_CONSUMING_OPS = {0, 1, 4, 7, 8}  # M, I, S, =, X


@dataclass
class CigarOperation:
    """Represents a single CIGAR operation with genomic coordinates."""
    op_code: int
    op_char: str
    length: int
    ref_start: int
    ref_end: int
    query_start: int
    query_end: int


@dataclass
class Deletion:
    """Represents a deletion detected from CIGAR."""
    ref_start: int
    ref_end: int
    size: int
    query_position: int  # Position in read where deletion occurs


@dataclass
class Insertion:
    """Represents an insertion detected from CIGAR."""
    ref_position: int  # Reference position before insertion
    size: int
    query_start: int
    query_end: int
    inserted_seq: Optional[str] = None


def parse_cigar_to_operations(read: pysam.AlignedSegment) -> List[CigarOperation]:
    """
    Parse CIGAR string into detailed operations with coordinates.

    Args:
        read: pysam AlignedSegment object

    Returns:
        List of CigarOperation objects with genomic coordinates
    """
    if read.cigartuples is None:
        return []

    operations = []
    ref_pos = read.reference_start
    query_pos = 0

    for op_code, length in read.cigartuples:
        op_char = CIGAR_OPS.get(op_code, '?')

        # Calculate end positions
        ref_end = ref_pos + length if op_code in REF_CONSUMING_OPS else ref_pos
        query_end = query_pos + length if op_code in QUERY_CONSUMING_OPS else query_pos

        operations.append(CigarOperation(
            op_code=op_code,
            op_char=op_char,
            length=length,
            ref_start=ref_pos,
            ref_end=ref_end,
            query_start=query_pos,
            query_end=query_end,
        ))

        # Update positions
        if op_code in REF_CONSUMING_OPS:
            ref_pos += length
        if op_code in QUERY_CONSUMING_OPS:
            query_pos += length

    return operations


def get_deletions_from_cigar(
    read: pysam.AlignedSegment,
    min_size: int = 1,
    target_region: Optional[Tuple[int, int]] = None
) -> List[Deletion]:
    """
    Extract deletions from CIGAR string.

    Args:
        read: pysam AlignedSegment object
        min_size: Minimum deletion size to report
        target_region: Optional (start, end) tuple to filter deletions

    Returns:
        List of Deletion objects
    """
    deletions = []
    operations = parse_cigar_to_operations(read)

    for op in operations:
        if op.op_code == 2 and op.length >= min_size:  # D = deletion
            # Check if deletion overlaps target region
            if target_region is not None:
                region_start, region_end = target_region
                if op.ref_end < region_start or op.ref_start > region_end:
                    continue

            deletions.append(Deletion(
                ref_start=op.ref_start,
                ref_end=op.ref_end,
                size=op.length,
                query_position=op.query_start,
            ))

    return deletions


def get_insertions_from_cigar(
    read: pysam.AlignedSegment,
    min_size: int = 1,
    target_region: Optional[Tuple[int, int]] = None
) -> List[Insertion]:
    """
    Extract insertions from CIGAR string.

    Args:
        read: pysam AlignedSegment object
        min_size: Minimum insertion size to report
        target_region: Optional (start, end) tuple to filter insertions

    Returns:
        List of Insertion objects
    """
    insertions = []
    operations = parse_cigar_to_operations(read)

    for op in operations:
        if op.op_code == 1 and op.length >= min_size:  # I = insertion
            # Check if insertion is within target region
            if target_region is not None:
                region_start, region_end = target_region
                if op.ref_start < region_start or op.ref_start > region_end:
                    continue

            # Get inserted sequence if available
            inserted_seq = None
            if read.query_sequence is not None:
                inserted_seq = read.query_sequence[op.query_start:op.query_end]

            insertions.append(Insertion(
                ref_position=op.ref_start,
                size=op.length,
                query_start=op.query_start,
                query_end=op.query_end,
                inserted_seq=inserted_seq,
            ))

    return insertions


def get_aligned_pairs(
    read: pysam.AlignedSegment,
    region_start: Optional[int] = None,
    region_end: Optional[int] = None
) -> List[Tuple[Optional[int], Optional[int], Optional[str]]]:
    """
    Get aligned pairs (query_pos, ref_pos, ref_base) for a read.

    Optionally filter to a specific reference region.

    Args:
        read: pysam AlignedSegment object
        region_start: Optional start of region to extract
        region_end: Optional end of region to extract

    Returns:
        List of (query_pos, ref_pos, ref_base) tuples
    """
    pairs = read.get_aligned_pairs(with_seq=True)

    if region_start is None and region_end is None:
        return pairs

    filtered_pairs = []
    for query_pos, ref_pos, ref_base in pairs:
        if ref_pos is None:
            continue
        if region_start is not None and ref_pos < region_start:
            continue
        if region_end is not None and ref_pos > region_end:
            continue
        filtered_pairs.append((query_pos, ref_pos, ref_base))

    return filtered_pairs


def count_mismatches_in_region(
    read: pysam.AlignedSegment,
    region_start: int,
    region_end: int
) -> Tuple[int, int, List[Tuple[int, str, str]]]:
    """
    Count matches and mismatches in a specific region.

    Args:
        read: pysam AlignedSegment object
        region_start: Start of region
        region_end: End of region

    Returns:
        Tuple of (n_matches, n_mismatches, mismatch_details)
        mismatch_details is list of (ref_pos, ref_base, query_base)
    """
    pairs = get_aligned_pairs(read, region_start, region_end)
    query_seq = read.query_sequence

    n_matches = 0
    n_mismatches = 0
    mismatch_details = []

    for query_pos, ref_pos, ref_base in pairs:
        if query_pos is None or ref_pos is None or ref_base is None:
            continue

        query_base = query_seq[query_pos] if query_seq else None

        if query_base and ref_base:
            if query_base.upper() == ref_base.upper():
                n_matches += 1
            else:
                n_mismatches += 1
                mismatch_details.append((ref_pos, ref_base.upper(), query_base.upper()))

    return n_matches, n_mismatches, mismatch_details


def get_mismatch_signature(read: pysam.AlignedSegment) -> str:
    """
    Generate a mismatch signature string for a read.

    Format: "pos1:ref>alt,pos2:ref>alt,..."
    Used for identifying reads with identical mismatch patterns.

    Args:
        read: pysam AlignedSegment object

    Returns:
        Mismatch signature string
    """
    _, _, mismatches = count_mismatches_in_region(
        read,
        read.reference_start,
        read.reference_end
    )

    if not mismatches:
        return ""

    sig_parts = [f"{pos}:{ref}>{alt}" for pos, ref, alt in sorted(mismatches)]
    return ",".join(sig_parts)


def get_soft_clip_count(read: pysam.AlignedSegment) -> int:
    """Get total soft-clipped bases from a read."""
    if read.cigartuples is None:
        return 0

    return sum(length for op, length in read.cigartuples if op == 4)


def get_total_indel_size(read: pysam.AlignedSegment) -> Tuple[int, int]:
    """Get total insertion and deletion sizes from a read.

    Returns:
        Tuple of (total_insertions, total_deletions) in bp
    """
    if read.cigartuples is None:
        return 0, 0

    insertions = sum(length for op, length in read.cigartuples if op == 1)
    deletions = sum(length for op, length in read.cigartuples if op == 2)

    return insertions, deletions
