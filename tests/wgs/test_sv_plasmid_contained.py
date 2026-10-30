"""Half B keeps donor-contained fragments only when the plasmids are known to be purged (planning/328).

Kevin R. Roy, 2026-10-30: "count contained read pairs when the backbone is empty ... a per clone/per dataset
option for trace". These tests pin the option (``keep_contained``), the plasmid-load measurement that decides
it per sample (``plasmid_load``), and the purge test (``plasmids_purged``).
"""

import math

import pysam
import pytest

from trace_crispr.wgs.sv_plasmid import (
    PlasmidLoad,
    SamRecord,
    classify_donor_reads,
    plasmid_load,
    plasmids_purged,
    poisson_upper,
)

CHROM = "chrI"
DSTART, DEND = 100, 200
BL = {"chrREDI": [(1, 1000)]}


def _rec(qname, flag, pos, cigar, rnext, pnext, tlen, mapq=60):
    return SamRecord(qname, flag, CHROM, pos, cigar, rnext, pnext, tlen, mapq=mapq)


def _contained():
    # both mates inside the widened donor [94, 206]: fragment 110..189
    return [_rec("c", 99, 110, "60M", "=", 130, 80), _rec("c", 147, 130, "60M", "=", 110, -80)]


def _classify(records, **kw):
    return classify_donor_reads(records, CHROM, DSTART, DEND, BL, dup_mode="standard", **kw)


def test_default_still_discards_a_contained_fragment():
    res = _classify(_contained())
    assert {k.qname for k in res.discarded} == {"c"}
    assert res.frags_plasmid == 1 and res.frags_contained_kept == 0


def test_keep_contained_retains_it_apart_from_genome():
    res = _classify(_contained(), keep_contained=True)
    assert {k.qname for k in res.included} == {"c"}
    assert res.discarded == set()
    assert res.frags_contained_kept == 1
    assert res.frags_genome == 0 and res.frags_plasmid == 0
    assert res.variant_cov == 0


def test_keep_contained_still_discards_a_donor_edge_clip():
    # contained fragment whose read is soft-clipped at the donor start: a junction, not a clean molecule
    recs = [_rec("j", 99, 100, "10S50M", "=", 120, 80), _rec("j", 147, 120, "60M", "=", 100, -80)]
    res = _classify(recs, keep_contained=True)
    assert {k.qname for k in res.discarded} == {"j"}
    assert res.frags_contained_kept == 0


def test_keep_contained_still_discards_a_cassette_mate():
    res = _classify([_rec("m", 97, 120, "50M", "chrREDI", 500, 0)], keep_contained=True)
    assert {k.qname for k in res.discarded} == {"m"}
    assert res.frags_contained_kept == 0


def test_keep_contained_does_not_touch_a_genomic_fragment():
    recs = [_rec("g", 99, 50, "75M", "=", 160, 185), _rec("g", 147, 160, "75M", "=", 50, -185)]
    res = _classify(recs, keep_contained=True)
    assert res.frags_genome == 1 and res.frags_contained_kept == 0


def test_parity_refuses_keep_contained():
    with pytest.raises(ValueError):
        classify_donor_reads(_contained(), CHROM, DSTART, DEND, BL, dup_mode="parity", keep_contained=True)


@pytest.mark.parametrize("n,expected", [(0, 3.689), (1, 5.572), (2, 7.225), (10, 18.390)])
def test_poisson_upper_matches_the_exact_limits(n, expected):
    assert math.isclose(poisson_upper(n), expected, abs_tol=0.002)


def _write_bam(path, nuclear_reads, plasmid_starts, donor_pair=False):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": "chrI", "LN": 100000}, {"SN": "pGD", "LN": 2000}]}
    reads = []

    def seg(name, ref, pos, flag=0, mate=None):
        a = pysam.AlignedSegment()
        a.query_name, a.flag, a.reference_id, a.reference_start = name, flag, ref, pos
        a.mapping_quality, a.cigarstring = 60, "50M"
        a.query_sequence, a.query_qualities = "A" * 50, pysam.qualitystring_to_array("I" * 50)
        if mate is not None:
            a.next_reference_id, a.next_reference_start = mate
        return a

    for i in range(nuclear_reads):
        reads.append(seg(f"n{i}", 0, 1000 + (i * 37) % 90000))
    for i, p in enumerate(plasmid_starts):
        reads.append(seg(f"p{i}", 1, p))
    if donor_pair:
        reads.append(seg("d", 0, 120, flag=1 | 0x40, mate=(1, 1500)))
    reads.sort(key=lambda a: (a.reference_id, a.reference_start))
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for a in reads:
            out.write(a)
    pysam.index(str(path))
    return str(path)


def test_plasmid_load_counts_backbone_outside_the_mask(tmp_path):
    bam = _write_bam(tmp_path / "s.bam", 2000, [10, 300, 900, 950, 1200])
    load = plasmid_load(bam, ["pGD"], mask=[("pGD", 851, 1050)], nuclear_contigs=["chrI"])
    assert load.plasmid_reads == 3            # 900 and 950 sit in the mask
    assert load.unmasked_length == 1800
    assert load.nuclear_reads == 2000
    assert math.isclose(load.copies, 3 / (2000 / 100000 * 1800))
    assert load.donor_backbone_pairs == 0


def test_purge_needs_depth_and_no_donor_backbone_pair(tmp_path):
    deep = plasmid_load(_write_bam(tmp_path / "deep.bam", 60000, []), ["pGD"], nuclear_contigs=["chrI"])
    shallow = plasmid_load(_write_bam(tmp_path / "shallow.bam", 500, []), ["pGD"], nuclear_contigs=["chrI"])
    paired = plasmid_load(_write_bam(tmp_path / "pair.bam", 60000, [], donor_pair=True), ["pGD"],
                          donor_windows=[("chrI", DSTART, DEND)], nuclear_contigs=["chrI"])
    assert deep.plasmid_reads == 0 and plasmids_purged(deep)          # 3.689 / 1200 expected reads per copy
    assert shallow.plasmid_reads == 0 and not plasmids_purged(shallow)  # zero reads at 10 per copy proves nothing
    assert paired.donor_backbone_pairs == 1 and not plasmids_purged(paired)


def test_plasmid_load_refuses_an_undeclared_contig(tmp_path):
    bam = _write_bam(tmp_path / "s.bam", 10, [])
    with pytest.raises(ValueError):
        plasmid_load(bam, ["pMissing"], nuclear_contigs=["chrI"])


def test_purged_rejects_a_nonpositive_threshold():
    load = PlasmidLoad(0, 1000, 100, 1000, 0, 100.0, 0.0, 0.0369)
    with pytest.raises(ValueError):
        plasmids_purged(load, max_copies=0)
