"""Unit tests for trace_crispr.wgs.donor_outcome (synthetic fixtures, no external data)."""

from pathlib import Path

import pytest

from trace_crispr.wgs.donor_outcome import (
    DESIGNED,
    SPONTANEOUS,
    SYNTHETIC_ERROR,
    DonorSpec,
    apply_variant,
    bucket_variant,
    levenshtein,
    orient_to_genome,
    precall_sample,
    reverse_complement,
)

# A synthetic donor window placed at chrT:11-20.
WT = "AAAAATTTTT"
DESIGN = "AACAATTTTT"  # intended SNV: donor-local pos 2 (genomic 13) A->C
STEP2 = "AACAATTGTT"  # design + a synthesis error at donor-local pos 7 (T->G)
REF = "N" * 10 + WT + "N" * 10  # chrT reference; donor at 1-based 11..20
DSTART, DEND = 11, 20


def _spec(sample="s1", wt=WT, design=DESIGN, step2=STEP2):
    return DonorSpec(sample, "chrT", DSTART, DEND, wt, design, step2)


# --------------------------- primitives ----------------------------------- #


def test_levenshtein_matches_known():
    assert levenshtein("AAAA", "AAAA") == 0
    assert levenshtein("AAAA", "AABA") == 1
    assert levenshtein("AAAA", "AAA") == 1  # deletion
    assert levenshtein("AAA", "AAAA") == 1  # insertion
    assert levenshtein("", "ACGT") == 4


def test_reverse_complement():
    assert reverse_complement("AAAAATTTTT") == "AAAAATTTTT"
    assert reverse_complement("ACGT") == "ACGT"
    assert reverse_complement("AACAATTTTT") == reverse_complement("AACAATTTTT")
    assert reverse_complement("ATGC") == "GCAT"
    with pytest.raises(ValueError):
        reverse_complement("AXGT")


def test_apply_variant_snv_ins_del():
    assert apply_variant("AAAAA", 2, "A", "C") == "AACAA"  # SNV
    assert apply_variant("AAAAA", 2, "A", "AG") == "AAAGAA"  # insertion
    assert apply_variant("AACAA", 2, "CA", "C") == "AACA"  # deletion


# --------------------------- bucketing ------------------------------------ #


def test_bucket_designed():
    donor_tmp = apply_variant(WT, 2, "A", "C")  # -> DESIGN exactly
    bucket, e2, e1 = bucket_variant(donor_tmp, WT, DESIGN, STEP2)
    assert bucket == DESIGNED
    assert (e2, e1) == (1, 0)  # wt is 1 from design; tmp is 0


def test_bucket_synthetic_error():
    # A variant that only moves toward the step2 (synthetic) donor, not the design.
    donor_tmp = apply_variant(WT, 7, "T", "G")  # WT pos7 T->G ; appears in STEP2 only
    bucket, _, _ = bucket_variant(donor_tmp, WT, DESIGN, STEP2)
    # e1 = lev(tmp, design): tmp="AAAAATTGTT" vs design "AACAATTTTT" -> 2 (pos2,pos7)
    # e2 = lev(wt, design) = 1 ; e1>=e2 ; e3 = lev(tmp, step2)=1 < e4 = lev(wt, step2)=2
    assert bucket == SYNTHETIC_ERROR


def test_bucket_spontaneous():
    # A variant unrelated to design or step2 (moves toward neither).
    donor_tmp = apply_variant(WT, 9, "T", "A")  # last base T->A, not in design/step2
    bucket, _, _ = bucket_variant(donor_tmp, WT, DESIGN, STEP2)
    assert bucket == SPONTANEOUS


def test_bucket_tie_is_not_designed():
    # If applying the variant leaves distance-to-design unchanged (e1==e2), strict
    # '<' means it is NOT designed.
    wt, design, step2 = "AAAA", "AAAA", "AAAA"
    donor_tmp = apply_variant(wt, 0, "A", "C")  # e1=1, e2=0 -> e1>=e2
    bucket, e2, e1 = bucket_variant(donor_tmp, wt, design, step2)
    assert bucket != DESIGNED
    assert e1 >= e2


# --------------------------- strand handling ------------------------------ #


def test_orient_to_genome_plus_strand_unchanged():
    spec = _spec()
    out = orient_to_genome(spec, REF)
    assert out.donor_wt == WT  # already matches reference window


def test_orient_to_genome_minus_strand_flips():
    # Use an asymmetric donor (NOT a reverse-complement palindrome) so the flip
    # is actually exercised.
    wt, design, step2 = "ACGTAACCGG", "ACCTAACCGG", "ACCTAACCGT"
    ref = "N" * 10 + wt + "N" * 10
    spec = DonorSpec(
        "s",
        "chrT",
        DSTART,
        DEND,
        reverse_complement(wt),
        reverse_complement(design),
        reverse_complement(step2),
    )
    out = orient_to_genome(spec, ref)
    assert out.donor_wt == wt
    assert out.donor_design == design
    assert out.donor_step2 == step2


def test_orient_to_genome_mismatch_raises():
    spec = _spec(wt="GGGGGGGGGG")
    with pytest.raises(ValueError):
        orient_to_genome(spec, REF)


# --------------------------- per-sample end-to-end ------------------------ #


def _write_vcf(path: Path, records):
    """records: list of (chrom, pos, ref, alt)."""
    with open(path, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, pos, ref, alt in records:
            fh.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t50\tPASS\t.\n")


def _write_gvcf(path: Path, records):
    """records: list of (chrom, pos, ref, alt_field, ref_depth)."""
    with open(path, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n")
        for chrom, pos, ref, alt, rd in records:
            fh.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT:AD:DP\t0:{rd},0:9\n")


def test_precall_installed(tmp_path):
    vcf = tmp_path / "s.fil.vcf"
    gvcf = tmp_path / "s.g.vcf"
    _write_vcf(vcf, [("chrT", 13, "A", "C")])  # the one designed SNV
    _write_gvcf(gvcf, [])
    out = precall_sample(_spec(), vcf, gvcf, REF)
    assert out.editing_outcome == "designed_variant_installed"
    assert out.designed_variants == ["chrT_13_A_to_C"]
    assert out.synthetic_errors == [] and out.spontaneous_variants == []


def test_precall_unedited_when_window_covered(tmp_path):
    vcf = tmp_path / "s.fil.vcf"
    gvcf = tmp_path / "s.g.vcf"
    _write_vcf(vcf, [])  # no variants
    # quantification window for 11..20 is [11,20]; cover all 10 with ref depth >=1
    _write_gvcf(gvcf, [("chrT", p, "A", "<NON_REF>", 8) for p in range(11, 21)])
    out = precall_sample(_spec(), vcf, gvcf, REF)
    assert out.editing_outcome == "designed_variant_unedited"


def test_precall_unknown_when_uncovered(tmp_path):
    vcf = tmp_path / "s.fil.vcf"
    gvcf = tmp_path / "s.g.vcf"
    _write_vcf(vcf, [])
    # cover only 5 of 10 window positions -> not unedited, not installed -> unknown
    _write_gvcf(gvcf, [("chrT", p, "A", "<NON_REF>", 8) for p in range(11, 16)])
    out = precall_sample(_spec(), vcf, gvcf, REF)
    assert out.editing_outcome == "unknown"


def test_precall_partial(tmp_path):
    # design needs the pos-13 SNV; supply a different in-window variant that moves
    # toward design's neighborhood but does not complete it -> partial OR unknown.
    # Use a 2-SNV design so one designed SNV alone is "partial".
    wt = "AAAAATTTTT"
    design = "ACAAATTTAT"  # two SNVs: pos1 A->C and pos8 T->A
    step2 = design
    spec = DonorSpec("s", "chrT", DSTART, DEND, wt, design, step2)
    vcf = tmp_path / "s.fil.vcf"
    gvcf = tmp_path / "s.g.vcf"
    _write_vcf(vcf, [("chrT", 12, "A", "C")])  # only the first designed SNV (pos1)
    _write_gvcf(gvcf, [])
    out = precall_sample(spec, vcf, gvcf, "N" * 10 + wt + "N" * 10)
    assert out.editing_outcome == "designed_variant_partially_installed"
    assert out.designed_variants == ["chrT_12_A_to_C"]
