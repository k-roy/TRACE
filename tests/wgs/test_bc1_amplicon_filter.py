"""Shengdi's amplicon telltale, and its Nextera read-through exemption at the read's 3' edge.

The telltale drops a read whose aligned start (end) lies in the left (right) amplicon
window with a soft clip or indel in its first 5 (last 7) query bases: the inner-primer
randomer. On untrimmed WGS input a genuine Tn5 read that runs off its fragment end is
soft-clipped too, into the Nextera mosaic end. That clip is at the read's 3' end, so it is
exempt there; the 5' edge, where the randomer sits, is never exempt.
"""
import pysam
import pytest

from trace_crispr.wgs.bc1_consensus import TrackConfig, is_amplicon_read, revcomp

MOSAIC_END = "CTGTCTCTTATACACATCT"
MOSAIC_TAIL = MOSAIC_END + "CCGAGCCCACGAGAC"  # read-through continues into the i7 adapter
TRUSEQ_TAIL = "AGATCGGAAGAGCGTCGTGTAG"
RANDOMER = "ACTTCA"
GENOMIC = "A" * 60
LEFT, RIGHT = (91, 102), (193, 219)
HEADER = pysam.AlignmentHeader.from_references(["chrII"], [1000])


def track(left=LEFT, right=RIGHT):
    return TrackConfig(name="genome", contig="chrII", bc_start=140, bc_end=159,
                       left_junction="ACCCTGTTAT", right_junction="GAGCATGACC", revcomp=True,
                       amplicon_left_window=left, amplicon_right_window=right)


def read(seq_ref, cigar, start1, reverse):
    """A primary proper-pair read; `seq_ref` is in reference orientation, as in a BAM."""
    r = pysam.AlignedSegment(HEADER)
    r.query_name = "read"
    r.query_sequence = seq_ref
    r.query_qualities = pysam.qualitystring_to_array("I" * len(seq_ref))
    r.flag = 83 if reverse else 99
    r.reference_id = r.next_reference_id = 0
    r.reference_start = start1 - 1
    r.next_reference_start = 300
    r.mapping_quality = 60
    r.cigarstring = cigar
    return r


def left_clipped(clip_ref, reverse):
    """Leading clip in reference orientation, aligned start 101 (inside LEFT)."""
    return read(clip_ref + GENOMIC, f"{len(clip_ref)}S{len(GENOMIC)}M", 101, reverse)


def right_clipped(clip_ref, reverse):
    """Trailing clip in reference orientation, aligned end 200 (inside RIGHT)."""
    return read(GENOMIC + clip_ref, f"{len(GENOMIC)}M{len(clip_ref)}S", 141, reverse)


def test_randomer_at_a_forward_reads_5prime_edge_is_an_amplicon():
    assert is_amplicon_read(left_clipped(RANDOMER, reverse=False), track())


def test_randomer_at_a_reverse_reads_5prime_edge_is_an_amplicon():
    # A reverse read's 5' end is its right edge; its randomer is revcomp'd in the BAM.
    assert is_amplicon_read(right_clipped(revcomp(RANDOMER), reverse=True), track())


def test_mosaic_end_read_through_at_a_reverse_reads_3prime_edge_is_exempt():
    # Sequencing orientation: genomic, then the mosaic end. In the BAM that tail is the
    # reverse complement, on the left, as a leading clip at the aligned start.
    r = left_clipped(revcomp(MOSAIC_TAIL), reverse=True)
    assert not is_amplicon_read(r, track())


def test_mosaic_end_read_through_at_a_forward_reads_3prime_edge_is_exempt():
    assert not is_amplicon_read(right_clipped(MOSAIC_TAIL, reverse=False), track())


@pytest.mark.parametrize("builder,clip", [
    (left_clipped, revcomp(MOSAIC_TAIL)),   # forward read: the left edge is its 5' end
    (right_clipped, MOSAIC_TAIL),            # reverse read: the right edge is its 5' end
])
def test_the_5prime_edge_is_never_exempt(builder, clip):
    reverse = builder is right_clipped
    assert is_amplicon_read(builder(clip, reverse=reverse), track())


def test_truseq_read_through_at_the_3prime_edge_is_still_an_amplicon():
    # An amplicon read running off its insert reads the partner randomer, then TruSeq.
    tail = revcomp(RANDOMER) + TRUSEQ_TAIL
    assert is_amplicon_read(right_clipped(tail, reverse=False), track())
    assert is_amplicon_read(left_clipped(revcomp(tail), reverse=True), track())


@pytest.mark.parametrize("aligned", [1, 3, 7])
def test_mosaic_end_partly_aligned_by_chance_is_exempt(aligned):
    # bwa aligns adapter bases that happen to match the reference, so the clip starts
    # mid-adapter; the aligned adapter bases sit just before it.
    seq = GENOMIC[:-aligned] + MOSAIC_END[:aligned] + MOSAIC_TAIL[aligned:]
    r = read(seq, f"{len(GENOMIC)}M{len(MOSAIC_TAIL) - aligned}S", 141, reverse=False)
    assert not is_amplicon_read(r, track())


def _mutate(s, positions):
    swap = {"A": "C", "C": "G", "G": "T", "T": "A"}
    return "".join(swap[b] if i in positions else b for i, b in enumerate(s))


@pytest.mark.parametrize("positions,exempt", [((), True), ((4,), True), ((4, 9), True),
                                              ((1, 4, 9, 12), False)])
def test_mismatch_tolerance_on_a_15_base_read_through(positions, exempt):
    clip = _mutate(MOSAIC_END[:15], positions)
    r = right_clipped(clip, reverse=False)
    assert is_amplicon_read(r, track()) is (not exempt)


@pytest.mark.parametrize("tail,exempt", [("CTGT", True), ("CAGT", False)])
def test_short_read_through_with_an_indel_and_no_clip(tail, exempt):
    # A 4-base run into the mosaic end that bwa aligned across a 1-base deletion: no clip,
    # but a defect in the last 7 bases at an aligned end inside RIGHT.
    r = read(GENOMIC + tail, f"{len(GENOMIC)}M1D{len(tail)}M", 136, reverse=False)
    assert r.reference_end in range(RIGHT[0], RIGHT[1] + 1)
    assert is_amplicon_read(r, track()) is (not exempt)


def test_a_read_with_both_a_5prime_randomer_and_3prime_read_through_is_an_amplicon():
    # The exemption is scoped to the 3' edge, so a randomer at the 5' edge still counts.
    seq = RANDOMER + GENOMIC + MOSAIC_TAIL
    r = read(seq, f"{len(RANDOMER)}S{len(GENOMIC)}M{len(MOSAIC_TAIL)}S", 101, reverse=False)
    assert is_amplicon_read(r, track())


def test_disabled_without_windows():
    t = track(left=None, right=None)
    assert not is_amplicon_read(left_clipped(RANDOMER, reverse=False), t)
    assert not is_amplicon_read(right_clipped(revcomp(RANDOMER), reverse=True), t)
