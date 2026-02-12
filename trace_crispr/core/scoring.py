"""
Alignment scoring and best alignment selection.

Implements triple-aligner consensus approach using BWA-MEM, BBMap, and minimap2.

Author: Kevin R. Roy
"""

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple
import pysam

from .cigar import get_soft_clip_count, get_total_indel_size


@dataclass
class AlignmentScore:
    """Unified alignment score for comparing aligners."""
    aligner: str
    ref_start: int = 0
    ref_end: int = 0
    mismatches: int = 0
    soft_clips: int = 0
    insertions: int = 0
    deletions: int = 0
    is_aligned: bool = False
    cigar: str = ""
    aligned_seq: str = ""
    md_tag: str = ""
    mate_unmapped: bool = False  # For paired-end: is the mate unmapped?

    @property
    def penalty_score(self) -> float:
        """Score for ranking alignments. Lower is better.
        Soft clips are penalized equally to mismatches."""
        if not self.is_aligned:
            return float('inf')
        return self.mismatches + self.soft_clips

    @property
    def priority(self) -> int:
        """Priority for tie-breaking: bwa=0 (highest), bbmap=1, minimap2=2"""
        priority_map = {'bwa': 0, 'bbmap': 1, 'minimap2': 2}
        return priority_map.get(self.aligner, 99)

    def __lt__(self, other: 'AlignmentScore') -> bool:
        """Compare alignments: lower penalty wins, then higher priority."""
        if self.penalty_score != other.penalty_score:
            return self.penalty_score < other.penalty_score
        return self.priority < other.priority


def score_alignment(read: pysam.AlignedSegment, aligner: str) -> AlignmentScore:
    """
    Score a single alignment from any aligner.

    Args:
        read: pysam AlignedSegment object
        aligner: Name of the aligner ('bwa', 'bbmap', 'minimap2')

    Returns:
        AlignmentScore object
    """
    if read.is_unmapped:
        return AlignmentScore(aligner=aligner, is_aligned=False)

    # Get soft clip count
    soft_clips = get_soft_clip_count(read)

    # Get indel sizes
    insertions, deletions = get_total_indel_size(read)

    # Get mismatches from MD tag or NM tag
    mismatches = 0
    if read.has_tag('NM'):
        # NM includes mismatches + insertions + deletions
        nm = read.get_tag('NM')
        mismatches = nm - insertions - deletions
    elif read.has_tag('MD'):
        # Parse MD tag for mismatches
        md_tag = read.get_tag('MD')
        mismatches = count_mismatches_from_md(md_tag)

    # Get CIGAR string
    cigar = read.cigarstring if read.cigarstring else ""

    # Get MD tag
    md_tag = read.get_tag('MD') if read.has_tag('MD') else ""

    return AlignmentScore(
        aligner=aligner,
        ref_start=read.reference_start,
        ref_end=read.reference_end,
        mismatches=mismatches,
        soft_clips=soft_clips,
        insertions=insertions,
        deletions=deletions,
        is_aligned=True,
        cigar=cigar,
        aligned_seq=read.query_sequence if read.query_sequence else "",
        md_tag=md_tag,
        mate_unmapped=read.mate_is_unmapped if read.is_paired else False,
    )


def count_mismatches_from_md(md_tag: str) -> int:
    """
    Count mismatches from an MD tag.

    MD tag format: [0-9]+([A-Z]|\\^[A-Z]+)[0-9]+...
    Numbers indicate matches, letters indicate mismatches,
    ^[A-Z]+ indicates deletions.
    """
    import re

    if not md_tag:
        return 0

    # Remove deletion markers (^[A-Z]+)
    md_clean = re.sub(r'\^[A-Z]+', '', md_tag)

    # Count remaining letters (mismatches)
    mismatches = sum(1 for c in md_clean if c.isalpha())

    return mismatches


def select_best_alignment(
    alignments: Dict[str, AlignmentScore]
) -> Tuple[str, AlignmentScore]:
    """
    Select the best alignment from multiple aligners.

    Selection criteria:
    1. Lower penalty score (mismatches + soft_clips)
    2. Ties go to higher priority aligner (bwa > bbmap > minimap2)

    Args:
        alignments: Dict mapping aligner name to AlignmentScore

    Returns:
        Tuple of (best_aligner_name, best_alignment_score)
    """
    if not alignments:
        return "", AlignmentScore(aligner="none", is_aligned=False)

    # Filter to aligned only
    aligned = {k: v for k, v in alignments.items() if v.is_aligned}

    if not aligned:
        # Return first unaligned
        first_key = next(iter(alignments))
        return first_key, alignments[first_key]

    # Find best
    best_aligner = min(aligned.keys(), key=lambda k: aligned[k])
    return best_aligner, aligned[best_aligner]


def select_best_alignment_paired(
    r1_alignments: Dict[str, AlignmentScore],
    r2_alignments: Dict[str, AlignmentScore]
) -> Tuple[str, AlignmentScore, AlignmentScore]:
    """
    Select the best alignment for paired-end reads.

    For paired-end, we sum the penalty scores from R1 and R2
    for each aligner and select the aligner with lowest combined score.

    Args:
        r1_alignments: Dict mapping aligner name to R1 AlignmentScore
        r2_alignments: Dict mapping aligner name to R2 AlignmentScore

    Returns:
        Tuple of (best_aligner_name, r1_score, r2_score)
    """
    # Get common aligners
    common_aligners = set(r1_alignments.keys()) & set(r2_alignments.keys())

    if not common_aligners:
        return "", AlignmentScore(aligner="none"), AlignmentScore(aligner="none")

    # Calculate combined scores
    combined_scores = {}
    for aligner in common_aligners:
        r1_score = r1_alignments[aligner]
        r2_score = r2_alignments[aligner]

        # Both must be aligned for a valid pair
        if not r1_score.is_aligned or not r2_score.is_aligned:
            combined_scores[aligner] = (float('inf'), 99)  # (penalty, priority)
        else:
            combined_penalty = r1_score.penalty_score + r2_score.penalty_score
            priority = min(r1_score.priority, r2_score.priority)
            combined_scores[aligner] = (combined_penalty, priority)

    # Find best aligner
    best_aligner = min(common_aligners, key=lambda k: combined_scores[k])

    return best_aligner, r1_alignments[best_aligner], r2_alignments[best_aligner]


@dataclass
class DeduplicationSignature:
    """Signature for deduplication of paired-end reads."""
    r1_ref_start: int
    r1_ref_end: int
    r2_ref_start: int
    r2_ref_end: int
    r1_cigar: str
    r2_cigar: str

    def __hash__(self):
        return hash((
            self.r1_ref_start, self.r1_ref_end,
            self.r2_ref_start, self.r2_ref_end,
            self.r1_cigar, self.r2_cigar
        ))

    def __eq__(self, other):
        if not isinstance(other, DeduplicationSignature):
            return False
        return (
            self.r1_ref_start == other.r1_ref_start and
            self.r1_ref_end == other.r1_ref_end and
            self.r2_ref_start == other.r2_ref_start and
            self.r2_ref_end == other.r2_ref_end and
            self.r1_cigar == other.r1_cigar and
            self.r2_cigar == other.r2_cigar
        )


def get_dedup_signature(
    r1_score: AlignmentScore,
    r2_score: AlignmentScore
) -> Optional[DeduplicationSignature]:
    """
    Create a deduplication signature for a paired-end read.

    Args:
        r1_score: R1 alignment score
        r2_score: R2 alignment score

    Returns:
        DeduplicationSignature or None if reads are not properly aligned
    """
    if not r1_score.is_aligned or not r2_score.is_aligned:
        return None

    return DeduplicationSignature(
        r1_ref_start=r1_score.ref_start,
        r1_ref_end=r1_score.ref_end,
        r2_ref_start=r2_score.ref_start,
        r2_ref_end=r2_score.ref_end,
        r1_cigar=r1_score.cigar,
        r2_cigar=r2_score.cigar,
    )
