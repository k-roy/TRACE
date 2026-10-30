"""Haplotype-level calls with the target codon judged apart from its synonymous helpers (planning/328, v7).

Kevin R. Roy, 2026-10-30: 2.2.2 is "an incomplete edit when a designed base has a wt position", with a separate
correct-target category when only synonymous helpers are wild type; amino-acid designs are scored on the
translation frame; "Syn designs can be treated liked missense for the syn helper mutations, which are distint
from the target codon"; "separate category, helpers only".
"""

import math

import pysam
import pytest

from trace_crispr.wgs import target_calls as tc
from trace_crispr.wgs.reconcile import IntendedAllele

PAD = "AAAAAAAAAA"
# M A K Q L F *  on the plus strand, CDS 11..31
CDS_SEQ = "ATG" "GCT" "AAA" "CAG" "CTG" "TTC" "TAA"
COMP = str.maketrans("ACGT", "TGCA")


@pytest.fixture
def fasta(tmp_path):
    plus = PAD + CDS_SEQ + PAD
    minus = (PAD + CDS_SEQ + PAD).translate(COMP)[::-1]
    path = tmp_path / "ref.fa"
    path.write_text(">chrP\n" + plus + "\n>chrM\n" + minus + "\n>chrBad\n" + PAD + "ATGTAGTAA" + PAD + "\n")
    pysam.faidx(str(path))
    with pysam.FastaFile(str(path)) as fa:
        yield fa


def _plus_missense(fa):
    # helpers: GCT->GCC (16 T>C), CTG->CTA (25 G>A); target CAG (Q) -> GCG (A): 20 C>G, 21 A>C
    alleles = [IntendedAllele(16, "T", "C"), IntendedAllele(20, "C", "G"), IntendedAllele(21, "A", "C"),
               IntendedAllele(25, "G", "A")]
    spec = tc.TargetSpec(tc.CDS("chrP", 11, 31, "+"), 20, 22)
    return tc.build_target_model(fa, alleles, spec)


def _plus_synonymous(fa):
    # target AAA (K) -> AAG (K): 19 A>G; helpers 16 T>C and 25 G>A
    alleles = [IntendedAllele(16, "T", "C"), IntendedAllele(19, "A", "G"), IntendedAllele(25, "G", "A")]
    return tc.build_target_model(fa, alleles, tc.TargetSpec(tc.CDS("chrP", 11, 31, "+"), 17, 19))


def call(model, frags):
    return tc.call_haplotypes(frags, model)


def test_model_separates_target_and_helpers(fasta):
    m = _plus_missense(fasta)
    assert m.target_wt_codon == "CAG" and m.target_design_codon == "GCG" and not m.synonymous_target
    assert m.helpers == frozenset({0, 3})
    assert _plus_synonymous(fasta).synonymous_target


@pytest.mark.parametrize("frags,nt,target", [
    (["DDDD"], "CE", "CE"),
    (["WDDW"], "IE", "CE_target"),          # target as designed, both helpers wild type
    (["DWWD"], "IE", "helpers_only"),       # helpers in, target wild type
    (["DWWD", "DWWD", "WWWW"], "IE", "helpers_only"),
    (["WWWW", "WWWW"], "NE", "NE"),
    (["WWWW"], "UN", "UN"),
    (["D..."], "CE", "CE_unconfirmed_target"),
    (["WDWW"], "IE", "IE"),                 # half the target codon: GAG (E), neither Q nor A
    (["DDDD"] + ["WWWW"] * 9, "MIX", "MIX"),
    (["DDDD"] * 9 + ["WWWW"], "CE", "CE"),  # one wild-type read pair among nine is what contamination gives
    (["WWWW"] * 4 + ["WOWW"], "NE", "NE"),  # an unexpected base is not a partial edit
    (["DDDD", "DWWD"], "IE", "MIX"),        # target-design and helpers-only molecules in one clone
])
def test_missense_calls(fasta, frags, nt, target):
    r = call(_plus_missense(fasta), frags)
    assert (r.call_nt, r.call_target) == (nt, target)


@pytest.mark.parametrize("frags,target", [(["DDD"], "CE"), (["WDW"], "CE_target"), (["DWD"], "helpers_only"),
                                          (["WWW", "WWW"], "NE")])
def test_synonymous_target_is_judged_on_its_codon(fasta, frags, target):
    assert call(_plus_synonymous(fasta), frags).call_target == target


def test_minus_strand_mirrors_plus(fasta):
    # the same design on the reverse-complemented contig: genomic positions mirror (p -> 42 - p), bases complement
    L = len(PAD + CDS_SEQ + PAD)
    def mirror(p, ref, alt):
        return IntendedAllele(L + 1 - p, ref.translate(COMP), alt.translate(COMP))
    alleles = sorted([mirror(16, "T", "C"), mirror(20, "C", "G"), mirror(21, "A", "C"), mirror(25, "G", "A")],
                     key=lambda a: a.pos)
    spec = tc.TargetSpec(tc.CDS("chrM", L + 1 - 31, L + 1 - 11, "-"), L + 1 - 22, L + 1 - 20)
    m = tc.build_target_model(fasta, alleles, spec)
    assert m.target_wt_codon == "CAG" and m.target_design_codon == "GCG"
    # alleles now run 25, 21, 20, 16 in genomic order: helper, target, target, helper
    assert call(m, ["WDDW"]).call_target == "CE_target"
    assert call(m, ["DWWD"]).call_target == "helpers_only"


def test_refuses_a_cds_that_does_not_translate(fasta):
    with pytest.raises(ValueError):
        tc.build_target_model(fasta, [IntendedAllele(14, "T", "A")], tc.TargetSpec(tc.CDS("chrBad", 11, 19, "+"), 14, 16))


def test_refuses_an_off_frame_target_codon(fasta):
    with pytest.raises(ValueError):
        tc.build_target_model(fasta, [IntendedAllele(20, "C", "G")], tc.TargetSpec(tc.CDS("chrP", 11, 31, "+"), 19, 21))


def test_binomial_tail():
    assert math.isclose(tc.binom_tail(1, 2, 0.054), 1 - 0.946 ** 2)
    assert tc.significant_minority(1, 9, 0.054, 0.01) and not tc.significant_minority(9, 1, 0.054, 0.01)


def test_fragment_haplotypes_map_observations():
    assert tc.fragment_haplotypes([["alt", "ref", None, "conflict", "T"]]) == ["DW..O"]
