"""Unit + golden-regression tests for :mod:`trace_crispr.wgs.bc1_consensus`.

The golden test is built from the real MAGESTIC-QTL WGS well 20250331_Tn5_plate_1D
sample_13 (genome track): a single genuine barcode read whose barcode sits in a
leading soft-clip, plus the on-contig bc1-amplicon proper-pair contaminants and the
cross-contig read-2 hops that surround it. It locks in, together, the three
empirically-earned behaviours — canonical orientation (genome barcode is read out as
the rev-comp of the truth), the soft-clip-extended overlap window (a right-anchored
read starting AT bc_end is kept), and the amplicon-telltale filter (proper-pair
on-contig amplicons are dropped; the genuine read is not). Ground truth from the
orthogonal REDI amplicon (SHA345/yT182): canonical bc1 == ``ATAAAAATGTTGCGATATCT``.
"""

import pysam

from trace_crispr.wgs.bc1_consensus import (
    BRIDGE_AMBIGUOUS,
    BRIDGE_EXACT,
    BRIDGE_HD1,
    BRIDGE_NONE,
    BridgeIndex,
    TrackConfig,
    _boundary_defect,
    call_sample,
    caller_confidence,
    extract_barcode,
    extract_track,
    hd1_neighbors,
    is_amplicon_read,
    is_read2_hop,
    ref_to_query_offset,
    revcomp,
    top_barcode,
)

# --------------------------------------------------------------------------- #
# Track + read fixtures                                                        #
# --------------------------------------------------------------------------- #

CONTIG = "yT177_after_bc1_integration"
HOP_CONTIG = "chrIV"
BC_START, BC_END = 10270, 10289  # 1-based inclusive (matches input.yaml)
LEFT_JUNC = "A" * 150
RIGHT_JUNC = "C" * 150

# Real sample_13 genome-track reads (FLAG, 1-based POS, MAPQ, CIGAR, mate-contig, SEQ).
GENUINE = (147, 10290, 40, "42S33M", "=",
           "AAGTTATGCACGTATCCTAGGGAGATATCGCAACATTTTTATGAGCATGACCTGTCGACGTCGTAGGAGCTCATC")
AMPLICON_LEFT = (99, 10234, 0, "4S37M34S", "=",
                 "TTCCGCATACATTATACGAAGTTATGCACGTATCCTAGGGGCTCCGTAGTTCCATTGGTGGCGCCTGACCTGTCG")
AMPLICON_RIGHT = (147, 10289, 0, "38S31M6S", "=",
                  "TTATGCACGTATCCTAGGGGCTCCGTAGTTCCATTGGTGGAGCATGACCTGTCGACGTCGTAGGAGCTCGTAAAA")
HOP_A = (161, 10233, 0, "3S37M35S", "chrIV",
         "ACGAGCATACATTATACGAAGTTATGCACGTATCCTAGGGTAGTAATAGCCAGTTTTTTCGAGCATGACCTGTCG")
HOP_B = (161, 10234, 1, "4S38M33S", "chrIV",
         "GCATGCATACATTATACGAAGTTATGCACGTATCCTAGGGGGAATTTTGAGGTACCATGAGAGCATGACCTGTCG")

TRUTH = "ATAAAAATGTTGCGATATCT"  # canonical bc1 (orthogonal REDI truth)
GENOME_FWD = "AGATATCGCAACATTTTTAT"  # == revcomp(TRUTH); how the genome contig reads it


def genome_track(amplicon: bool = True) -> TrackConfig:
    return TrackConfig(
        name="genome", contig=CONTIG, bc_start=BC_START, bc_end=BC_END,
        left_junction=LEFT_JUNC, right_junction=RIGHT_JUNC, revcomp=True,
        amplicon_left_window=(10224, 10235) if amplicon else None,
        amplicon_right_window=(10318, 10326) if amplicon else None,
    )


def _header() -> pysam.AlignmentHeader:
    return pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": CONTIG, "LN": 14710}, {"SN": HOP_CONTIG, "LN": 1000}],
    })


def _segment(read, header):
    flag, pos1, mapq, cig, mate, seq = read
    a = pysam.AlignedSegment(header)
    a.query_name = f"r{pos1}_{flag}"
    a.flag = flag
    a.reference_id = 0
    a.reference_start = pos1 - 1
    a.mapping_quality = mapq
    a.cigarstring = cig
    a.query_sequence = seq
    a.next_reference_id = 0 if mate == "=" else 1
    a.next_reference_start = pos1 - 1 if mate == "=" else 100
    return a


# --------------------------------------------------------------------------- #
# Pure helpers                                                                 #
# --------------------------------------------------------------------------- #


def test_revcomp_preserves_n():
    assert revcomp("AAAACCCCGGGGTTTTACGT") == "ACGTAAAACCCCGGGGTTTT"
    assert revcomp("ACGTN") == "NACGT"


def test_to_canonical_only_revcomps_when_flagged():
    g = TrackConfig("g", CONTIG, BC_START, BC_END, LEFT_JUNC, RIGHT_JUNC, revcomp=True)
    p = TrackConfig("p", CONTIG, BC_START, BC_END, LEFT_JUNC, RIGHT_JUNC, revcomp=False)
    assert g.to_canonical(GENOME_FWD) == TRUTH
    assert p.to_canonical(TRUTH) == TRUTH


def test_boundary_defect_first_and_last():
    h = _header()
    # leading 4S -> defect within first 5 query bases
    assert _boundary_defect(_segment(AMPLICON_LEFT, h), 5, from_end=False) is True
    # trailing 6S -> defect within last 7
    assert _boundary_defect(_segment(AMPLICON_RIGHT, h), 7, from_end=True) is True
    # genuine read: clean 33M at the right boundary -> no last-7 defect
    assert _boundary_defect(_segment(GENUINE, h), 7, from_end=True) is False


# --------------------------------------------------------------------------- #
# ref_to_query_offset + extract_barcode                                        #
# --------------------------------------------------------------------------- #


def test_ref_to_query_offset_into_leading_softclip():
    """Right-anchored genuine read: bc_end maps to the match start; the barcode is
    the 20 query bases just before it, recovered out of the leading soft-clip."""
    h = _header()
    seg = _segment(GENUINE, h)
    q_end = ref_to_query_offset(seg, BC_END)  # 0-based one-past-last barcode base
    assert q_end == 42  # 42S leading soft-clip ends here
    assert seg.query_sequence[q_end - 20:q_end] == GENOME_FWD


def test_extract_barcode_right_anchor_genuine_read():
    seg = _segment(GENUINE, _header())
    ext = extract_barcode(seg, genome_track())
    assert ext is not None
    assert ext.anchor == "right"
    assert ext.query_barcode == GENOME_FWD
    assert ext.barcode == TRUTH  # canonicalised via revcomp


def test_extract_barcode_left_anchor_synthetic():
    """Left-anchored read: 40M ends exactly at the flank/barcode boundary, barcode in
    the trailing soft-clip (the common genome-track geometry)."""
    h = _header()
    bc_fwd = "TTTTGGGGCCCCAAAACGTA"
    query = "A" * 40 + bc_fwd + "C" * 15  # 75bp: flank(40) + barcode(20) + right(15)
    a = pysam.AlignedSegment(h)
    a.query_name = "left_synth"
    a.flag = 99
    a.reference_id = 0
    a.reference_start = BC_START - 1 - 40  # 40M ends at bc_start-1 (0-based)
    a.mapping_quality = 30
    a.cigarstring = "40M35S"
    a.query_sequence = query
    a.next_reference_id = 0
    a.next_reference_start = a.reference_start
    ext = extract_barcode(a, genome_track())
    assert ext is not None
    assert ext.anchor == "left"
    assert ext.query_barcode == bc_fwd
    assert ext.barcode == revcomp(bc_fwd)
    assert ext.flank_mismatches == 0  # aligned 40 'A' match LEFT_JUNC


# --------------------------------------------------------------------------- #
# Contaminant predicates                                                       #
# --------------------------------------------------------------------------- #


def test_is_amplicon_read_drops_proper_pair_amplicons_keeps_genuine():
    h, t = _header(), genome_track()
    assert is_amplicon_read(_segment(AMPLICON_LEFT, h), t) is True   # start window + first-5 defect
    assert is_amplicon_read(_segment(AMPLICON_RIGHT, h), t) is True  # end window + last-7 defect
    # genuine read ends INSIDE the right window (10322) but has no last-7 defect -> kept
    assert is_amplicon_read(_segment(GENUINE, h), t) is False


def test_is_amplicon_read_disabled_without_windows():
    h = _header()
    assert is_amplicon_read(_segment(AMPLICON_LEFT, h), genome_track(amplicon=False)) is False


def test_is_read2_hop_on_mate_off_contig():
    h, t = _header(), genome_track()
    assert is_read2_hop(_segment(HOP_A, h), t) is True
    assert is_read2_hop(_segment(GENUINE, h), t) is False


# --------------------------------------------------------------------------- #
# Golden end-to-end: extract_track over the sample_13 slice                    #
# --------------------------------------------------------------------------- #


def _write_bam(tmp_path, reads):
    h = _header()
    unsorted = tmp_path / "u.bam"
    with pysam.AlignmentFile(str(unsorted), "wb", header=h) as out:
        for r in reads:
            out.write(_segment(r, h))
    sorted_bam = tmp_path / "s.bam"
    pysam.sort("-o", str(sorted_bam), str(unsorted))
    pysam.index(str(sorted_bam))
    return str(sorted_bam)


def test_extract_track_golden_recovers_truth_and_drops_contaminants(tmp_path):
    reads = [HOP_A, AMPLICON_LEFT, HOP_B, AMPLICON_RIGHT, GENUINE]
    bam = _write_bam(tmp_path, reads)
    with pysam.AlignmentFile(bam, "rb") as af:
        call = extract_track(af, genome_track(amplicon=True))
    support = call.barcode_support()
    # the genuine read survives; the two on-contig amplicons + two hops are dropped
    assert set(support) == {TRUTH}
    assert support[TRUTH] == 1
    assert call.n_amplicon == 2
    assert call.n_contaminant >= 4  # 2 amplicon + 2 hops


def test_extract_track_without_amplicon_filter_lets_contaminant_win(tmp_path):
    """Documents WHY the amplicon filter is needed: with it off, the on-contig
    amplicon proper-pair (support 2) out-votes the genuine read (support 1)."""
    reads = [HOP_A, AMPLICON_LEFT, HOP_B, AMPLICON_RIGHT, GENUINE]
    bam = _write_bam(tmp_path, reads)
    with pysam.AlignmentFile(bam, "rb") as af:
        call = extract_track(af, genome_track(amplicon=False))
    support = call.barcode_support()
    assert call.n_amplicon == 0
    assert TRUTH in support  # genuine still extracted...
    # ...but a foreign barcode from the amplicon pair has >= its support
    assert max(support.values()) >= support[TRUTH]
    assert len(support) >= 2


def test_window_shift_breaks_the_call(tmp_path):
    """+-1 bp frame shift must stop recovering the truth (the empirical shift control)."""
    reads = [GENUINE]
    bam = _write_bam(tmp_path, reads)
    shifted = TrackConfig(
        name="genome", contig=CONTIG, bc_start=BC_START + 1, bc_end=BC_END + 1,
        left_junction=LEFT_JUNC, right_junction=RIGHT_JUNC, revcomp=True,
    )
    with pysam.AlignmentFile(bam, "rb") as af:
        call = extract_track(af, shifted)
    assert TRUTH not in call.barcode_support()


# --------------------------------------------------------------------------- #
# Bridge: exact / HD1 / ambiguous / none + in-set-HD1 flag                     #
# --------------------------------------------------------------------------- #

B_EXACT = "AAAACCCCGGGGTTTTACGT"
B_EXACT_NB = "CAAACCCCGGGGTTTTACGT"  # HD1 neighbour of B_EXACT
B_AMBIG = "TTTTGGGGCCCCAAAATGCA"
B_SOLO = "ACGTACGTACGTACGTAAAA"


def _bridge():
    return BridgeIndex.from_records([
        (B_EXACT, "oligoA", 100), (B_EXACT, "oligoA", 50),  # single oligo, support 150
        (B_EXACT_NB, "oligoNB", 7),                          # in-set HD1 neighbour of B_EXACT
        (B_AMBIG, "oligoX", 30), (B_AMBIG, "oligoY", 40),    # two oligos -> ambiguous
        (B_SOLO, "oligoZ", 20),                              # isolated -> HD1 target
    ])


def test_bridge_exact_with_inset_hd1_flag():
    r = _bridge().bridge(B_EXACT)
    assert r.tier == BRIDGE_EXACT
    assert r.oligo == "oligoA"
    assert r.support == 150
    assert r.has_inset_hd1_neighbor is True
    assert "INSET_HD1_NEIGHBOR" in r.flags


def test_bridge_ambiguous_when_multiple_oligos():
    r = _bridge().bridge(B_AMBIG)
    assert r.tier == BRIDGE_AMBIGUOUS
    assert r.n_oligos == 2


def test_bridge_hd1_unique():
    query = B_SOLO[:-1] + "C"  # HD1 from B_SOLO only
    r = _bridge().bridge(query)
    assert r.tier == BRIDGE_HD1
    assert r.oligo == "oligoZ"
    assert r.matched_barcode == B_SOLO


def test_bridge_none_when_no_hit_or_neighbor():
    r = _bridge().bridge("GGGGGGGGGGGGGGGGGGGG")
    assert r.tier == BRIDGE_NONE
    assert r.oligo is None


def test_bridge_hd1_disabled():
    query = B_SOLO[:-1] + "C"
    r = _bridge().bridge(query, allow_hd1=False)
    assert r.tier == BRIDGE_NONE


def test_hd1_neighbors_count():
    assert len(set(hd1_neighbors("ACGTACGTACGTACGTACGT"))) == 60


# --------------------------------------------------------------------------- #
# top_barcode + caller_confidence + call_sample end-to-end                     #
# --------------------------------------------------------------------------- #


def test_top_barcode_and_confidence(tmp_path):
    bam = _write_bam(tmp_path, [GENUINE])
    with pysam.AlignmentFile(bam, "rb") as af:
        call = extract_track(af, genome_track())
    top = top_barcode(call)
    assert top.barcode == TRUTH
    assert top.support == 1 and top.both_anchors is False
    assert caller_confidence(top) == "low"  # single-read single-anchor


def test_call_sample_end_to_end(tmp_path):
    reads = [HOP_A, AMPLICON_LEFT, AMPLICON_RIGHT, GENUINE]
    bam = _write_bam(tmp_path, reads)
    bridge = BridgeIndex.from_records([(TRUTH, "oligoNNS", 932)])
    with pysam.AlignmentFile(bam, "rb") as af:
        call = call_sample(af, "well_X", genome_track(), plasmid_track=None, bridge=bridge)
    assert call.sample == "well_X"
    assert call.track == "genome"
    assert call.barcode == TRUTH
    assert call.bridge_tier == BRIDGE_EXACT
    assert call.oligo == "oligoNNS"
    assert call.caller_confidence == "low"   # genuine read is a single right-anchored molecule
    assert call.confidence == "medium"       # low caller x EXACT bridge -> medium
