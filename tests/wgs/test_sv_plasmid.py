"""Unit tests for :mod:`trace_crispr.wgs.sv_plasmid` (synthetic fixtures, CI-able).

These exercise the fidelity-critical Half B primitives (the duplicate hack, the
M+D span, the interval predicates, the full branch matrix of
:func:`classify_donor_reads`), plus Half A/C and the reference builder. The
empirical bit-for-bit parity against Shengdi's reference ``*.target_filtering.anno``
lives in ``tests/wgs/parity_clean_donor_reads.py`` (needs Sherlock data + samtools).
"""

import pytest

from trace_crispr.wgs.sv_plasmid import (
    SamRecord,
    _is_contained,
    _is_mapped_to_cassette,
    _is_overlapped,
    _is_softclipped_near_junction,
    build_integrated_plasmid_contig,
    classify_donor_reads,
    decode_duplicate,
    integrated_plasmid_junction,
    one_read_inside_large_deletion,
    ref_span_md,
    rev_comp,
    within_read_large_deletions,
)

# --------------------------------------------------------------------------- #
# decode_duplicate — mechanical ground truth from the Perl if_duplicate        #
# --------------------------------------------------------------------------- #
#
# Captured by running Shengdi's exact ``if_duplicate`` over flags 0..4095:
#   - 1536 flags return true;
#   - flag 1024 (the *real* 0x400 duplicate bit) returns FALSE — the bug;
#   - in [1024, 2047] exactly the 512 ODD flags are true (substr reads bit 0);
#   - in [2048, 4095] exactly the 1024 flags with 0x2 set are true (bit 1).
_PERL_N_TRUE = 1536


def test_decode_duplicate_parity_matches_perl_over_all_4096_flags():
    trues = [f for f in range(4096) if decode_duplicate(f, mode="parity")]
    assert len(trues) == _PERL_N_TRUE
    # No flag below 1024 is ever flagged (binary string shorter than 11 chars).
    assert all(f >= 1024 for f in trues)
    # Per-range structure (independently derived from the Perl substr arithmetic).
    assert sum(1 for f in range(1024, 2048) if decode_duplicate(f, "parity")) == 512
    assert sum(1 for f in range(2048, 4096) if decode_duplicate(f, "parity")) == 1024
    assert all((f % 2 == 1) for f in range(1024, 2048) if decode_duplicate(f, "parity"))
    assert all((f & 0x2) for f in range(2048, 4096) if decode_duplicate(f, "parity"))


@pytest.mark.parametrize(
    "flag,expected",
    [
        (0, False),
        (1023, False),
        (1024, False),  # THE BUG: real duplicate flag reads as not-duplicate
        (1025, True),
        (1026, False),
        (1027, True),
        (2048, False),
        (2049, False),
        (2050, True),
        (2051, True),
        (4094, True),
        (4095, True),
    ],
)
def test_decode_duplicate_parity_spot_checks(flag, expected):
    assert decode_duplicate(flag, mode="parity") is expected


@pytest.mark.parametrize(
    "flag,expected",
    [(0, False), (1024, True), (1025, True), (3, False), (1024 | 0x2, True), (0x800, False)],
)
def test_decode_duplicate_standard_is_correct_0x400(flag, expected):
    assert decode_duplicate(flag, mode="standard") is expected


def test_decode_duplicate_rejects_unknown_mode():
    with pytest.raises(ValueError):
        decode_duplicate(0, mode="nope")


# --------------------------------------------------------------------------- #
# ref_span_md — M + D only                                                      #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    "cigar,span",
    [
        ("100M", 100),
        ("10S90M", 90),  # soft clip does not extend reference span
        ("5H20M5S", 20),  # hard + soft clip excluded
        ("50M2D48M", 100),  # deletion extends span
        ("30M5I30M", 60),  # insertion does NOT extend span
        ("*", 0),
    ],
)
def test_ref_span_md(cigar, span):
    assert ref_span_md(cigar) == span


def test_ref_span_md_rejects_skipped_region_like_perl():
    # The Perl regex does not recognise N (it would die); we surface that as ValueError.
    with pytest.raises(ValueError):
        ref_span_md("20M100N20M")


# --------------------------------------------------------------------------- #
# Interval predicates                                                           #
# --------------------------------------------------------------------------- #


def test_is_overlapped():
    assert _is_overlapped(100, 200, 150, 250)
    assert _is_overlapped(100, 200, 200, 300)  # touching counts as overlap
    assert not _is_overlapped(100, 200, 201, 300)
    assert not _is_overlapped(100, 200, 0, 99)


def test_is_contained():
    assert _is_contained(110, 190, 100, 200)
    assert _is_contained(100, 200, 100, 200)  # equal bounds contained
    assert not _is_contained(90, 190, 100, 200)
    assert not _is_contained(110, 210, 100, 200)


def test_is_softclipped_near_junction():
    # leading soft clip near the left donor boundary
    assert _is_softclipped_near_junction(100, 200, 100, "10S90M", 6)
    assert _is_softclipped_near_junction(100, 200, 104, "10S90M", 6)
    assert not _is_softclipped_near_junction(100, 200, 120, "10S90M", 6)
    # trailing soft clip near the right donor boundary: pos + span(M+D) vs right
    assert _is_softclipped_near_junction(100, 200, 110, "90M10S", 6)  # 110+90=200
    assert not _is_softclipped_near_junction(100, 200, 80, "90M10S", 6)  # 80+90=170
    # no soft clip
    assert not _is_softclipped_near_junction(100, 200, 100, "100M", 6)


def test_load_blacklist_permissive_parsing(tmp_path):
    from trace_crispr.wgs.sv_plasmid import load_blacklist

    f = tmp_path / "bl.txt"
    # tab rows with optional 4th NAME column, plus one inert space-delimited row
    f.write_text(
        "chrXV\t721500\t723000\tHIS3\n"
        "yT177_before_bc1_integration\t1\t15295\n"
        "yT177_after_bc1_integration 1   14710\n"  # space-delimited -> inert
    )
    bl = load_blacklist(str(f))
    assert bl["chrXV"] == [(721500, 723000)]  # NAME column ignored
    assert bl["yT177_before_bc1_integration"] == [(1, 15295)]
    # the mangled row is kept under the whole-string key with an inert (0,0) range
    assert bl["yT177_after_bc1_integration 1   14710"] == [(0, 0)]
    assert not _is_mapped_to_cassette("yT177_after_bc1_integration 1   14710", 5000, bl)


def test_is_mapped_to_cassette_strict_bounds():
    bl = {"chrREDI": [(100, 200)], "chrV": [(500, 600)]}
    assert _is_mapped_to_cassette("chrREDI", 150, bl)
    assert not _is_mapped_to_cassette("chrREDI", 100, bl)  # strict: not at start
    assert not _is_mapped_to_cassette("chrREDI", 200, bl)  # strict: not at end
    assert not _is_mapped_to_cassette("chrOTHER", 150, bl)


# --------------------------------------------------------------------------- #
# classify_donor_reads — full branch matrix                                     #
# --------------------------------------------------------------------------- #

CHROM = "chrI"
DSTART, DEND = 100, 200  # v_coord = 150.0
BL = {"chrREDI": [(1, 1000)]}


def _rec(qname, flag, rname, pos, cigar, rnext, pnext, tlen):
    return SamRecord(qname, flag, rname, pos, cigar, rnext, pnext, tlen)


def _classify(records, dup_mode="parity"):
    return classify_donor_reads(records, CHROM, DSTART, DEND, BL, dup_mode=dup_mode)


def test_non_overlapping_read_ignored():
    # read [300,350] disjoint from donor [100,200]
    res = _classify([_rec("r", 99, CHROM, 300, "50M", "=", 300, 50)])
    assert res.overlapped == set()


def test_wrong_chromosome_ignored():
    res = _classify([_rec("r", 99, "chrII", 120, "50M", "=", 120, 50)])
    assert res.overlapped == set()


def test_near_mate_contained_is_plasmid():
    # fragment [110,189] contained in widened donor [94,206] -> discard/plasmid
    res = _classify([_rec("r", 99, CHROM, 110, "80M", "=", 110, 80)])
    assert res.discarded == {"r"}
    assert res.unique_plasmid == {"r"}
    assert res.frags_plasmid == 1 and res.frags_genome == 0


def test_near_mate_softclip_near_junction_is_chimeric():
    # fragment not contained (tlen 250), leading soft clip at donor_start -> chimeric
    res = _classify([_rec("r", 99, CHROM, 100, "10S90M", "=", 100, 250)])
    assert res.included == {"r"}
    assert res.unique_chimeric == {"r"}
    assert res.redi_translocation == {"r"}
    # pos 100 < v_coord 150 < read_end 190 -> variant_cov increments
    assert res.variant_cov == 1


def test_near_mate_clean_is_genome_and_variant_cov():
    # [120,180] not contained (tlen 300), no soft clip, spans midpoint -> genome
    res = _classify([_rec("r", 99, CHROM, 120, "60M", "=", 120, 300)])
    assert res.unique_genome == {"r"}
    assert res.variant_cov == 1
    assert res.frags_genome == 1


def test_distal_blacklist_read_exceeds_donor_is_chimeric():
    # mate distal on chrREDI (blacklist), read [100,219] exceeds widened donor -> chimeric
    res = _classify([_rec("r", 99, CHROM, 100, "120M", "chrREDI", 500, 0)])
    assert res.unique_chimeric == {"r"}
    assert res.redi_translocation == {"r"}
    # distal path never increments variant_cov
    assert res.variant_cov == 0


def test_distal_blacklist_read_contained_is_plasmid():
    # mate distal on chrREDI, read [120,169] contained in widened donor -> plasmid
    res = _classify([_rec("r", 99, CHROM, 120, "50M", "chrREDI", 500, 0)])
    assert res.unique_plasmid == {"r"}
    assert res.discarded == {"r"}


def test_distal_non_blacklist_is_genome():
    # mate distal on a non-blacklist chrom -> genome + variant_cov
    res = _classify([_rec("r", 99, CHROM, 110, "60M", "chrOTHER", 500, 0)])
    assert res.unique_genome == {"r"}
    assert res.variant_cov == 1


def test_far_tlen_same_chrom_takes_distal_path():
    # rnext '=' but |tlen| >= 1000 -> distal path; mate coord not in blacklist -> genome
    res = _classify([_rec("r", 99, CHROM, 110, "60M", "=", 5000, 5000)])
    assert res.unique_genome == {"r"}


def test_unmapped_mate_is_discarded_only():
    res = _classify([_rec("r", 99, CHROM, 130, "40M", "*", 0, 0)])
    assert res.discarded == {"r"}
    assert res.unique_plasmid == set() and res.unique_genome == set()


def test_duplicate_excluded_from_counts_but_kept_in_sets():
    # flag 2050 (0x800|0x2) is a "duplicate" under the parity hack -> not counted,
    # but the read is still classified (overlapped/included).
    res = _classify([_rec("r", 2050, CHROM, 120, "60M", "=", 120, 300)])
    assert res.included == {"r"}
    assert res.unique_genome == set()
    assert res.variant_cov == 0


def test_duplicate_standard_mode_counts_non_0x400_read():
    # flag 2050 has NO 0x400 bit, so the correct test treats it as not-a-duplicate
    # and it IS counted -- the parity-vs-standard contrast on the same flag.
    res = classify_donor_reads(
        [_rec("r", 2050, CHROM, 120, "60M", "=", 120, 300)],
        CHROM,
        DSTART,
        DEND,
        BL,
        dup_mode="standard",
    )
    assert res.unique_genome == {"r"}
    assert res.variant_cov == 1


def test_cross_category_double_count():
    # a pair whose two reads land in different categories is counted in both columns
    read1 = _rec("p", 99, CHROM, 120, "60M", "chrOTHER", 500, 0)  # genome
    read2 = _rec("p", 147, CHROM, 100, "120M", "chrREDI", 500, 0)  # chimeric
    res = _classify([read1, read2])
    assert "p" in res.unique_genome
    assert "p" in res.unique_chimeric
    assert res.frags_genome == 1 and res.frags_chimeric == 1


def test_variant_cov_counts_reads_not_fragments():
    # two non-duplicate reads of the same pair both span the midpoint -> +2
    read1 = _rec("p", 99, CHROM, 120, "60M", "chrOTHER", 500, 0)
    read2 = _rec("p", 147, CHROM, 130, "40M", "chrOTHER", 500, 0)  # [130,170] spans 150
    res = _classify([read1, read2])
    assert res.variant_cov == 2
    assert res.frags_genome == 1  # but only one fragment (read-name) for the count


# --------------------------------------------------------------------------- #
# Half C: mate-geometry + within-read large deletions                          #
# --------------------------------------------------------------------------- #


def test_one_read_inside_large_deletion_nested_and_wide():
    # read2 [2000,2100] nested in read1 [1000,5000]; span 4001 > 101 + 1000
    assert one_read_inside_large_deletion(1000, 5000, 2000, 2100)
    # symmetric: read1 nested in read2
    assert one_read_inside_large_deletion(2000, 2100, 1000, 5000)


def test_one_read_inside_large_deletion_requires_nesting():
    # overlap but not nested
    assert not one_read_inside_large_deletion(1000, 2000, 1500, 2500)


def test_one_read_inside_large_deletion_requires_wide_gap():
    # nested but the pair span barely exceeds the inner read
    assert not one_read_inside_large_deletion(1000, 1500, 1100, 1200)


class _FakeRead:
    """Minimal stand-in for a pysam AlignedSegment for cigar.get_deletions_from_cigar."""

    def __init__(self, cigartuples, reference_start=1000):
        self.cigartuples = cigartuples
        self.reference_start = reference_start
        self.query_sequence = None


def test_within_read_large_deletions_strict_cutoff():
    read = _FakeRead([(0, 50), (2, 60), (0, 50)])  # 60 bp deletion
    dels = within_read_large_deletions(read, cutoff=50)
    assert len(dels) == 1 and dels[0].size == 60
    # strict: a 60 bp deletion is not "large" when cutoff == 60
    assert within_read_large_deletions(read, cutoff=60) == []


# --------------------------------------------------------------------------- #
# Reference construction                                                         #
# --------------------------------------------------------------------------- #


def test_rev_comp_iupac():
    assert rev_comp("ACGT") == "ACGT"
    assert rev_comp("AACG") == "CGTT"
    assert rev_comp("ACGTN") == "NACGT"
    # ambiguity codes round-trip (R<->Y, etc.)
    assert rev_comp(rev_comp("RYMKVBHD")) == "RYMKVBHD"


def test_build_integrated_plasmid_contig_plus_strand():
    extended = "AAAA" + "xx" + "TTTT"  # bp_to_extend=4 -> flanks AAAA / TTTT
    contig = build_integrated_plasmid_contig(extended, "CCC", "GG", "+", bp_to_extend=4)
    assert contig == "AAAA" + "CCC" + "GG" + "CCC" + "TTTT"


def test_build_integrated_plasmid_contig_minus_strand_revcomps_segment():
    extended = "AAAA" + "xx" + "TTTT"
    contig = build_integrated_plasmid_contig(extended, "CCC", "GG", "-", bp_to_extend=4)
    assert contig == "AAAA" + "CCC" + "CC" + "CCC" + "TTTT"  # rev_comp("GG") == "CC"


def test_integrated_plasmid_junction():
    # bp_to_extend 1000 + 129 bp donor == 1129 (the parity literal)
    assert integrated_plasmid_junction("N" * 129, bp_to_extend=1000) == 1129
