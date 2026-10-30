"""Half B keeps a dovetailed genomic fragment that runs past the donor (SGTC planning/328 §4g).

The molecule is the N89N clone's `...0:0230` pair, in its own coordinates: the forward mate is aligned from
174051 (11S64M, its clip carrying the design), the reverse mate from 174045 (75M, to 174119). The insert is
shorter than two reads, so the reverse mate starts LEFT of the forward mate, and bwa's TLEN of 69 does not span
the reverse read's alignment. The fragment ends at 174119, 10 bp past the donor end, in chromosome sequence a
plasmid donor template cannot carry, so it is genomic and must be kept.
"""
from trace_crispr.wgs import sv_plasmid as sv

DSTART, DEND = 173981, 174109


def pair(mc=True):
    fwd = sv.SamRecord("frag", 99, "chrXIV", 174051, "11S64M", "=", 174045, 69, 60, "rg",
                       "75M" if mc else "")
    rev = sv.SamRecord("frag", 147, "chrXIV", 174045, "75M", "=", 174051, -69, 60, "rg",
                       "11S64M" if mc else "")
    return [fwd, rev]


def test_dovetailed_genomic_fragment_past_the_donor_is_kept():
    res = sv.classify_donor_reads(pair(), "chrXIV", DSTART, DEND, {})
    key = sv.FragmentKey("rg", "frag")
    assert key in res.included and key not in res.discarded and key in res.unique_genome


def test_without_a_mate_cigar_the_tlen_rule_still_applies():
    res = sv.classify_donor_reads(pair(mc=False), "chrXIV", DSTART, DEND, {})
    assert sv.FragmentKey("rg", "frag") in res.discarded


def test_a_fragment_inside_the_donor_is_still_discarded_as_a_template():
    fwd = sv.SamRecord("t", 163, "chrXIV", 174051, "17S25M", "=", 174051, 25, 60, "rg", "17S25M")
    rev = sv.SamRecord("t", 83, "chrXIV", 174051, "17S25M", "=", 174051, -25, 60, "rg", "17S25M")
    res = sv.classify_donor_reads([fwd, rev], "chrXIV", DSTART, DEND, {})
    assert sv.FragmentKey("rg", "t") in res.discarded


def test_parse_sam_line_reads_the_mate_cigar():
    line = "q\t147\tchrXIV\t174045\t60\t75M\t=\t174051\t-69\tSEQ\tQUAL\tMC:Z:11S64M\tRG:Z:rg"
    assert sv.parse_sam_line(line).mate_cigar == "11S64M"
