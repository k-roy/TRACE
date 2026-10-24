"""Molecule expectations fixed independently of alignment through barcode N20."""
import json
import random
from pathlib import Path

import pysam
import pytest

from trace_crispr.wgs.bc1_consensus import (
    TrackConfig, call_sample, extract_barcode, extract_track, revcomp,
)

BC = "ACGTCAGTTCGATGACCTGA"
LEFT = "GGTCATGCTC"
RIGHTS = ("ATAACAGGGT", "CCCTAGGATA")


def record(seq, cigar, start=90, flag=99, name="molecule"):
    header = pysam.AlignmentHeader.from_references(["cassette"], [1000])
    read = pysam.AlignedSegment(header)
    read.query_name = name
    read.query_sequence = seq
    read.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
    read.flag = flag
    read.reference_id = read.next_reference_id = 0
    read.reference_start = start
    read.next_reference_start = 180
    read.template_length = 140
    read.mapping_quality = 60
    read.cigarstring = cigar
    return read


def bam_file(tmp_path, reads):
    path = tmp_path / "reads.bam"
    with pysam.AlignmentFile(str(path), "wb", header=reads[0].header) as out:
        for read in sorted(reads, key=lambda r: r.reference_start):
            out.write(read)
    pysam.index(str(path))
    return path


@pytest.mark.parametrize("right", RIGHTS)
@pytest.mark.parametrize("canonical_reverse", [False, True])
@pytest.mark.parametrize("flag", [99, 83])
@pytest.mark.parametrize("cigar", ["40M", "10M5I10M5D15M", "10M20I20D10M"])
def test_both_cassettes_both_strands_ignore_barcode_alignment(right, canonical_reverse, flag, cigar):
    left = LEFT
    if canonical_reverse:
        left, right = revcomp(right), revcomp(left)
    raw = revcomp(BC) if canonical_reverse else BC
    track = TrackConfig("genome", "cassette", 101, 120, left, right, canonical_reverse)
    result = extract_barcode(record(left + raw + right, cigar, flag=flag), track)
    assert result is not None
    assert result.barcode == BC
    assert result.query_barcode == raw


@pytest.mark.parametrize("right", RIGHTS)
@pytest.mark.parametrize("anchor,cigar,start", [
    ("left", "10M20S", 90), ("right", "20S10M", 120),
    ("left", "10M20S5H", 90), ("right", "5H20S10M", 120),
])
def test_indexed_bam_one_flank_and_hard_plus_soft_clips(tmp_path, right, anchor, cigar, start):
    track = TrackConfig("genome", "cassette", 101, 120, LEFT, right)
    seq = LEFT + BC if anchor == "left" else BC + right
    path = bam_file(tmp_path, [record(seq, cigar, start)])
    with pysam.AlignmentFile(str(path)) as bam:
        output = call_sample(bam, "well_A1", track)
    assert output.sample == "well_A1"
    assert output.barcode == BC
    assert output.support == 1


@pytest.mark.parametrize("anchor", ["left", "right"])
@pytest.mark.parametrize("available,expected", [(3, None), (4, BC), (9, BC)])
def test_partial_flank_requires_four_exact_adjacent_bases(anchor, available, expected):
    track = TrackConfig("genome", "cassette", 101, 120, LEFT, RIGHTS[0])
    if anchor == "left":
        read = record(LEFT[-available:] + BC, f"{available}M20S", 100 - available)
    else:
        read = record(BC + RIGHTS[0][:available], f"20S{available}M", 120)
    result = extract_barcode(read, track)
    assert (result.barcode if result else None) == expected


@pytest.mark.parametrize("width", [19, 21, 30])
def test_both_flanks_with_wrong_barcode_width_refuse(width):
    seq = LEFT + (BC + "ACGTACGTAC")[:width] + RIGHTS[0]
    # Align both real constant flanks to their configured reference positions.
    middle = f"{width}I20D"
    track = TrackConfig("genome", "cassette", 101, 120, LEFT, RIGHTS[0])
    assert extract_barcode(record(seq, "10M" + middle + "10M"), track) is None


def test_wrong_width_opposite_flank_in_softclip_also_refuses():
    track = TrackConfig("genome", "cassette", 101, 120, LEFT, RIGHTS[0])
    read = record(LEFT + BC + "ACGTACGTAC" + RIGHTS[0], "10M40S")
    assert extract_barcode(read, track) is None


def test_two_complete_possible_frames_are_ambiguous():
    track = TrackConfig("genome", "cassette", 101, 120, LEFT, RIGHTS[0])
    # A deletion/insertion in the constant flank yields two independently
    # alignment-supported starts; neither may be chosen by list order.
    seq = LEFT + BC + RIGHTS[0] + LEFT + BC + RIGHTS[0]
    read = record(seq, "5M40I5M20I20D10M")
    assert extract_barcode(read, track) is None


def test_motif_elsewhere_without_locus_anchor_cannot_call():
    track = TrackConfig("genome", "cassette", 101, 120, LEFT, RIGHTS[0])
    assert extract_barcode(record(LEFT + BC + RIGHTS[0], "40M", start=200), track) is None


def test_wrong_anchor_sequence_cannot_call():
    track = TrackConfig("genome", "cassette", 101, 120, LEFT, RIGHTS[0])
    assert extract_barcode(record("T" * 10 + BC, "10M20S"), track) is None


@pytest.mark.parametrize("anchor", ["left", "right"])
def test_opposite_motif_inside_barcode_is_not_a_second_flank(anchor):
    track = TrackConfig("genome", "cassette", 101, 120, LEFT, RIGHTS[0])
    barcode = "ACGTA" + (RIGHTS[0] if anchor == "left" else LEFT) + "TGCAT"
    read = record(LEFT + barcode, "10M20S") if anchor == "left" else record(barcode + RIGHTS[0], "20S10M", 120)
    output = extract_barcode(read, track)
    assert output is not None
    assert output.barcode == barcode


@pytest.mark.parametrize("width", [19, 21, 30])
@pytest.mark.parametrize("mutated_side", ["left", "right"])
def test_wrong_width_with_mutated_aligned_opposite_flank_refuses(width, mutated_side):
    track = TrackConfig("genome", "cassette", 101, 120, LEFT, RIGHTS[0])
    left = "TGTCATGCTC" if mutated_side == "left" else LEFT
    right = "TTAACAGGGT" if mutated_side == "right" else RIGHTS[0]
    seq = left + (BC + "ACGTACGTAC")[:width] + right
    assert extract_barcode(record(seq, f"10M{width}I20D10M"), track) is None


def test_one_read_with_two_flanks_is_not_two_anchors_or_molecules(tmp_path):
    track = TrackConfig("genome", "cassette", 101, 120, LEFT, RIGHTS[0])
    path = bam_file(tmp_path, [record(LEFT + BC + RIGHTS[0], "40M")])
    with pysam.AlignmentFile(str(path)) as bam:
        output = call_sample(bam, "well", track)
    assert output.support == 1
    assert output.both_anchors is False
    assert output.caller_confidence == "low"


@pytest.mark.parametrize("width,last,expected", [(19, "A", None), (19, "T", None), (20, "A", BC), (20, "T", None)])
def test_one_aligned_opposite_base_checks_geometry_but_cannot_anchor_alone(width, last, expected):
    track = TrackConfig("genome", "cassette", 101, 120, LEFT, RIGHTS[0])
    read = record(LEFT + BC[:width] + last, f"10M{width}I20D1M")
    result = extract_barcode(read, track)
    assert (result.barcode if result else None) == expected


def test_validated_anchor_labels_feed_existing_pair_confidence_rule(tmp_path):
    track = TrackConfig("genome", "cassette", 101, 120, LEFT, RIGHTS[0])
    first = record(LEFT + BC, "10M20S", name="pair")
    second = record(LEFT[-3:] + BC + RIGHTS[0], "33M", start=97, flag=147, name="pair")
    first.next_reference_start = 97
    second.next_reference_start = 90
    first.template_length = 100
    second.template_length = -100
    # Three left bases are below the explicit minimum; the right flank supplies
    # the actual validated anchor, even though reference_start is left of N20.
    assert extract_barcode(second, track).anchor == "right"
    path = bam_file(tmp_path, [first, second])
    with pysam.AlignmentFile(str(path)) as bam:
        output = call_sample(bam, "well", track)
    assert output.support == 1
    assert output.both_anchors is True
    assert output.caller_confidence == "high"


def test_hard_clipped_bases_cannot_be_recovered(tmp_path):
    track = TrackConfig("genome", "cassette", 101, 120, LEFT, RIGHTS[0])
    path = bam_file(tmp_path, [record(LEFT, "10M20H")])
    with pysam.AlignmentFile(str(path)) as bam:
        assert call_sample(bam, "empty", track).barcode == ""


def test_softclipped_constant_tail_can_be_validated():
    track = TrackConfig("genome", "cassette", 101, 120, LEFT, RIGHTS[0])
    assert extract_barcode(record(LEFT + BC, "5M25S"), track).barcode == BC


def test_configured_cassette_reference_must_match(tmp_path):
    path = tmp_path / "cassette.fa"
    path.write_text(">cassette\n" + "A" * 90 + LEFT + "N" * 20 + RIGHTS[1] + "A" * 50 + "\n")
    pysam.faidx(str(path))
    with pysam.FastaFile(str(path)) as fasta:
        TrackConfig("genome", "cassette", 101, 120, LEFT, RIGHTS[1]).validate_against_reference(fasta)
        with pytest.raises(ValueError, match="right junction"):
            TrackConfig("genome", "cassette", 101, 120, LEFT, RIGHTS[0]).validate_against_reference(fasta)


def test_500_seeded_variable_region_cigars_preserve_molecular_barcode():
    rng = random.Random(20260919)
    for case_number in range(500):
        expected = "".join(rng.choice("ACGT") for _ in range(20))
        left, right = LEFT, RIGHTS[case_number % 2]
        canonical_reverse = bool(case_number % 3)
        if canonical_reverse:
            left, right = revcomp(right), revcomp(left)
        raw = revcomp(expected) if canonical_reverse else expected
        q = r = 0
        operations = [(0, 10)]
        while q < 20 or r < 20:
            choices = ([0] if q < 20 and r < 20 else []) + ([1] if q < 20 else []) + ([2] if r < 20 else [])
            op = rng.choice(choices)
            q += op in (0, 1)
            r += op in (0, 2)
            if operations[-1][0] == op:
                operations[-1] = (op, operations[-1][1] + 1)
            else:
                operations.append((op, 1))
        if operations[-1][0] == 0:
            operations[-1] = (0, operations[-1][1] + 10)
        else:
            operations.append((0, 10))
        read = record(left + raw + right, "40M")
        read.cigartuples = operations
        track = TrackConfig("genome", "cassette", 101, 120, left, right, canonical_reverse)
        output = extract_barcode(read, track)
        assert output is not None, (case_number, read.cigarstring)
        assert output.barcode == expected, (case_number, read.cigarstring)


REAL = json.loads((Path(__file__).parent / "fixtures/bc1_real_reads.json").read_text())


@pytest.mark.parametrize("case", REAL, ids=[r["id"] for r in REAL])
def test_frozen_real_reads_through_indexed_bam(tmp_path, case):
    c = case
    expected = c["expected_canonical"]
    raw = revcomp(expected) if c["revcomp"] else expected
    assert c["left_flank"] + raw + c["right_flank"] in c["seq"]
    header = pysam.AlignmentHeader.from_references([c["contig"]], [c["bc_end1"] + 1000])
    read = pysam.AlignedSegment(header)
    read.query_name = c["qname"]
    read.query_sequence = c["seq"]
    read.query_qualities = c["qualities"]
    read.flag = c["flag"]
    read.reference_id = read.next_reference_id = 0
    read.reference_start = c["pos1"] - 1
    read.next_reference_start = c["mate_pos1"] - 1
    read.template_length = c["template_length"]
    read.cigarstring = c["cigar"]
    read.mapping_quality = c["mapq"]
    track = TrackConfig(c["track"], c["contig"], c["bc_start1"], c["bc_end1"],
                        c["left_flank"], c["right_flank"], c["revcomp"], min_mapq=0)
    # This frozen replay tests molecular boundaries, including two MAPQ<20
    # reads. Production locus-identity MAPQ refusal has separate tests.
    path = bam_file(tmp_path, [read])
    with pysam.AlignmentFile(str(path)) as bam:
        assert call_sample(bam, c["sample"], track).barcode == expected
