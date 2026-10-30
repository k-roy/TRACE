"""Local realignment at a localized edit site (SGTC planning/328).

Every molecule below is specified before the implementation is called. The reference is a fixed pseudo-random
sequence so a realignment is not degenerate. The design is six substitutions packed into 16 bp, dense enough
that a genome-wide aligner soft-clips a read across it, which is the case the aligned-base reader cannot see.
"""
import random

import pysam
import pytest

from trace_crispr.wgs import reconcile as r
from trace_crispr.wgs import local_realign as lr

random.seed(328)
REF = "".join(random.choice("ACGT") for _ in range(1200))
CHROM = "chrT"
# six designed substitutions at 1-based 601..616, every third base
DESIGN_POS = [601, 604, 607, 610, 613, 616]
SWAP = {"A": "C", "C": "G", "G": "T", "T": "A"}
ALLELES = [(p, REF[p - 1], SWAP[REF[p - 1]]) for p in DESIGN_POS]


def edited(seq_start, length, upto=None):
    """Reference sequence from 1-based `seq_start`, carrying the design (all alleles, or those < `upto`)."""
    s = list(REF[seq_start - 1: seq_start - 1 + length])
    for p, _, alt in ALLELES:
        if upto is not None and p >= upto:
            continue
        if seq_start <= p < seq_start + length:
            s[p - seq_start] = alt
    return "".join(s)


def bam_with(tmp_path, records):
    fa = tmp_path / "ref.fa"
    fa.write_text(f">{CHROM}\n{REF}\n")
    pysam.faidx(str(fa))
    header = pysam.AlignmentHeader.from_references([CHROM], [len(REF)])
    reads = []
    for d in records:
        a = pysam.AlignedSegment(header)
        a.query_name = d["name"]
        a.query_sequence = d["seq"]
        a.flag = d.get("flag", 0)
        a.reference_id = 0
        a.reference_start = d["start"] - 1
        a.cigarstring = d["cigar"]
        a.mapping_quality = d.get("mapq", 60)
        q = [35] * len(d["seq"])
        for i in d.get("lowq", ()):
            q[i] = 5
        a.query_qualities = pysam.qualitystring_to_array("".join(chr(x + 33) for x in q))
        reads.append(a)
    path = tmp_path / "t.bam"
    with pysam.AlignmentFile(str(path), "wb", header=header) as f:
        for a in sorted(reads, key=lambda x: x.reference_start):
            f.write(a)
    pysam.index(str(path))
    return pysam.AlignmentFile(str(path), "rb"), pysam.FastaFile(str(fa))


def site():
    return r.EditSite(CHROM, 571, 646, "oligo_t", tuple(r.IntendedAllele(*a) for a in ALLELES))


def counts(bam, fa, alleles=ALLELES, **kw):
    c, unmapped, _ = lr.realigned_allele_counts(bam, fa, CHROM, alleles, **kw)
    return [(x["alt"], x["ref"], x["other"]) for x in c], unmapped


def test_clipped_read_across_dense_design_is_read_by_realignment(tmp_path):
    # aligned 541-600, the edited 16 bp + 4 more soft-clipped (as bwa does on a dense design)
    seq = edited(541, 80)
    bam, fa = bam_with(tmp_path, [dict(name="clip", start=541, seq=seq, cigar="60M20S")])
    got, _ = counts(bam, fa)
    assert got == [(1, 0, 0)] * 6
    old = [r.intended_allele_evidence(bam, fa, CHROM, r.IntendedAllele(*a)) for a in ALLELES]
    assert [e.depth for e in old] == [0] * 6          # the aligned-base reader sees nothing
    call = r.score_edit_site(bam, fa, site())
    assert call.status == r.SITE_FULL and "LOCAL_REALIGNMENT" in call.flags
    assert r.score_edit_site(bam, fa, site(), evidence="aligned").status == r.SITE_NOCOV


def test_left_clip_and_reverse_strand_read(tmp_path):
    # aligned from 618 onward, its first 22 bases (597..617, carrying the design) soft-clipped; reverse strand
    seq = edited(597, 80)
    bam, fa = bam_with(tmp_path, [dict(name="lclip", start=618, seq=seq, cigar="21S59M", flag=16)])
    got, _ = counts(bam, fa)
    assert got == [(1, 0, 0)] * 6


def test_wild_type_read_is_ref(tmp_path):
    bam, fa = bam_with(tmp_path, [dict(name="wt", start=561, seq=REF[560:640], cigar="80M")])
    got, _ = counts(bam, fa)
    assert got == [(0, 1, 0)] * 6


def test_partial_incorporation_is_kept_per_allele(tmp_path):
    seq = edited(561, 80, upto=608)                       # alleles at 601, 604, 607 installed; 610+ wild type
    bam, fa = bam_with(tmp_path, [dict(name=f"p{i}", start=561, seq=seq, cigar="80M") for i in range(4)])
    got, _ = counts(bam, fa)
    assert got == [(4, 0, 0)] * 3 + [(0, 4, 0)] * 3
    assert r.score_edit_site(bam, fa, site()).status == r.SITE_INCOMPLETE


def test_backbone_chimera_clipped_at_donor_edge_is_not_counted(tmp_path):
    random.seed(7)
    backbone = "".join(random.choice("ACGT") for _ in range(40))
    # genomic 531..570 then plasmid backbone: the clip projects over the designed span but is not genome sequence
    bam, fa = bam_with(tmp_path, [dict(name="chim", start=531, seq=REF[530:570] + backbone, cigar="40M40S")])
    got, _ = counts(bam, fa)
    assert got == [(0, 0, 0)] * 6


def test_low_quality_designed_base_is_not_counted(tmp_path):
    seq = edited(561, 80)
    bam, fa = bam_with(tmp_path, [dict(name="lq", start=561, seq=seq, cigar="80M", lowq=(601 - 561,))])
    got, _ = counts(bam, fa)
    assert got[0] == (0, 0, 0) and got[1:] == [(1, 0, 0)] * 5


def test_deletion_inside_the_clip_is_realigned_not_projected(tmp_path):
    # the read carries the design AND a 2-bp deletion at 591-592, inside its clipped part; naive projection of
    # the clip onto the reference would put every designed base 2 bp out of frame
    left = REF[540:590]                                   # 541..590
    right = edited(593, 40)                               # 593..632 with the design
    bam, fa = bam_with(tmp_path, [dict(name="del", start=541, seq=left + right, cigar="50M40S")])
    got, _ = counts(bam, fa)
    assert got == [(1, 0, 0)] * 6


def test_fragment_mates_count_once_and_conflicts_drop(tmp_path):
    seq = edited(561, 80)
    recs = [dict(name="f", start=561, seq=seq, cigar="80M", flag=99), dict(name="f", start=565, seq=edited(565, 76), cigar="76M", flag=147),
            dict(name="g", start=561, seq=seq, cigar="80M", flag=99), dict(name="g", start=561, seq=REF[560:640], cigar="80M", flag=147)]
    bam, fa = bam_with(tmp_path, recs)
    got, _ = counts(bam, fa)
    assert got == [(1, 0, 0)] * 6                         # f counts once; g's mates disagree and drop


def test_length_changing_allele_by_haplotype_comparison(tmp_path):
    pos = 700
    dele = (pos, REF[pos - 1: pos + 2], REF[pos - 1])     # delete two bases after the anchor
    carrier = REF[649:pos] + REF[pos + 2: pos + 52]       # 650..700 then 703..752
    bam, fa = bam_with(tmp_path, [dict(name="d1", start=650, seq=carrier, cigar="51M2D50M"),
                                  dict(name="w1", start=650, seq=REF[649:750], cigar="101M")])
    got, _ = counts(bam, fa, alleles=[dele])
    assert got == [(1, 1, 0)]


def test_align_scores_like_bwa_and_clips():
    score, cols = lr.align("ACGTACGTAC", "TTACGTACGTACTT")
    assert score == 10 and len(cols) == 10
    s2, c2 = lr.align("ACGTACGTAC" + "GGGGG", "TTACGTACGTACTTTTTTT")
    assert s2 == 10 - lr.CLIP and len(c2) == 10           # a trailing mismatch run is clipped, not aligned


@pytest.mark.parametrize("flag", [0x100, 0x200, 0x400, 0x800])
def test_nonprimary_duplicate_qcfail_reads_supply_no_evidence(tmp_path, flag):
    bam, fa = bam_with(tmp_path, [dict(name="bad", start=561, seq=edited(561, 80), cigar="80M", flag=flag),
                                  dict(name="ok", start=561, seq=REF[560:640], cigar="80M")])
    got, _ = counts(bam, fa)
    assert got == [(0, 1, 0)] * 6


def test_mapq_floor_and_unknown_mapq_are_excluded_and_counted(tmp_path):
    bam, fa = bam_with(tmp_path, [dict(name="lo", start=561, seq=edited(561, 80), cigar="80M", mapq=5),
                                  dict(name="unk", start=561, seq=edited(561, 80), cigar="80M", mapq=255)])
    c, _, _ = lr.realigned_allele_counts(bam, fa, CHROM, ALLELES, min_mapq=20)
    assert all((x["depth"], x["low_mapq_reads"], x["unknown_mapq_reads"]) == (0, 1, 1) for x in c)
    c0, _, _ = lr.realigned_allele_counts(bam, fa, CHROM, ALLELES, min_mapq=0)
    assert all((x["alt"], x["unknown_mapq_reads"]) == (1, 1) for x in c0)


def test_placed_unmapped_mate_is_counted_not_used(tmp_path):
    bam, fa = bam_with(tmp_path, [dict(name="u", start=561, seq=edited(561, 80), cigar="80M", flag=0x4 | 0x1 | 0x40)])
    got, unmapped = counts(bam, fa)
    assert got == [(0, 0, 0)] * 6 and unmapped == 1


def test_homopolymer_placement_keeps_the_original_alignment():
    target = "G" * 10 + "A" * 30 + "G" * 10
    s1, cols = lr.align("A" * 10, target, expect_end=29)
    assert s1 == 10 and [t for _, t in cols] == list(range(20, 30))


def test_a_read_slid_onto_a_design_inside_a_repeat_is_not_realigned(tmp_path):
    # a 200-A run with the design A500C inside it; the read lies wholly in the run and carries a sequencing-error
    # C at 501. Sliding it one base left would put that C on the designed base at a perfect score, so the
    # realignment moves the read; the anchor guard must refuse that and fall back to the aligned reading (ref).
    ref = REF[:400] + "A" * 200 + REF[600:]
    fa = tmp_path / "hp.fa"
    fa.write_text(f">{CHROM}\n{ref}\n")
    pysam.faidx(str(fa))
    header = pysam.AlignmentHeader.from_references([CHROM], [len(ref)])
    a = pysam.AlignedSegment(header)
    seq = "A" * 30 + "C" + "A" * 49                       # 471..550, C at 501
    a.query_name, a.query_sequence, a.flag, a.reference_id, a.reference_start = "noise", seq, 0, 0, 470
    a.cigarstring, a.mapping_quality = "80M", 60
    a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
    path = tmp_path / "hp.bam"
    with pysam.AlignmentFile(str(path), "wb", header=header) as f:
        f.write(a)
    pysam.index(str(path))
    with pysam.AlignmentFile(str(path)) as bam, pysam.FastaFile(str(fa)) as fasta:
        evidence, tally = r.realigned_intended_evidence(bam, fasta, CHROM, (r.IntendedAllele(500, "A", "C"),))
    assert (evidence[0].alt_count, evidence[0].ref_count) == (0, 1)
    assert tally["reads_realign_displaced"] == 1 and tally["observations_from_fallback"] == 1
