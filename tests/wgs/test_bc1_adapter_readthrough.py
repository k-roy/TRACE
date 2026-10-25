"""A barcode slice that runs into Nextera read-through is not a barcode.

On untrimmed WGS input, a Tn5 fragment that ends inside the N20 reads on into the mosaic end,
and the anchored 20-base slice then carries adapter bases. On plate 5031 this produced calls
such as `CTGTCTCTTATACACATCTC` (all adapter) and `AGAGACAG` + 12 genuine bases (15 counted
reads in 11 wells). A genuine barcode is followed by its constant flank, never by the mosaic
end, so the refusal keys on a mosaic end that starts inside the slice and continues past it.
"""
import pysam

from trace_crispr.wgs.bc1_consensus import TrackConfig, extract_barcode, revcomp

BC = "ACGTCAGTTCGATGACCTGA"
LEFT = "GGTCATGCTC"
RIGHT = "ATAACAGGGT"
MOSAIC_END = "CTGTCTCTTATACACATCT"
I7_TAIL = "CCGAGCCCACGAGAC"   # read 1 continues into the i7 adapter
I5_TAIL = "GACGCTGCCGACGA"    # read 2 continues into the i5 adapter


def record(seq, cigar, start, flag):
    header = pysam.AlignmentHeader.from_references(["cassette"], [1000])
    read = pysam.AlignedSegment(header)
    read.query_name = "molecule"
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


def track():
    return TrackConfig("genome", "cassette", 101, 120, LEFT, RIGHT)


def test_forward_read_whose_fragment_ends_inside_the_barcode_is_refused():
    seq = LEFT + BC[:12] + MOSAIC_END + I7_TAIL
    assert extract_barcode(record(seq, f"10M{len(seq) - 10}S", 90, 99), track()) is None


def test_reverse_read_whose_slice_is_entirely_adapter_is_refused():
    # Sequencing orientation: genomic (right flank onward), then the mosaic end at the
    # barcode's right edge. In the BAM that adapter is reverse-complemented on the left.
    clip = revcomp(MOSAIC_END + I5_TAIL)
    seq = clip + RIGHT
    assert extract_barcode(record(seq, f"{len(clip)}S10M", 120, 83), track()) is None


def test_read_through_with_two_sequencing_errors_is_still_refused():
    noisy = MOSAIC_END[:4] + "A" + MOSAIC_END[5:12] + "G" + MOSAIC_END[13:]
    seq = LEFT + BC[:12] + noisy + I7_TAIL
    assert extract_barcode(record(seq, f"10M{len(seq) - 10}S", 90, 99), track()) is None


def test_genuine_barcode_followed_by_its_flank_is_kept():
    seq = LEFT + BC + RIGHT
    result = extract_barcode(record(seq, "40M", 90, 99), track())
    assert result is not None and result.barcode == BC


def test_genuine_barcode_containing_a_mosaic_end_like_10mer_is_kept():
    # The slice holds CTGTCTCTTA, but what follows is the constant flank, not the mosaic end.
    lookalike = "GACTGTCTCTTAGGCATAAC"
    seq = LEFT + lookalike + RIGHT
    result = extract_barcode(record(seq, "40M", 90, 99), track())
    assert result is not None and result.barcode == lookalike


def test_mosaic_end_starting_one_base_before_the_slice_is_refused():
    # Plate 5031 s35: the fragment ended one base into the right flank, and bwa aligned the
    # mosaic end's first base (C) to the flank's first base (G on the other strand). The
    # adapter therefore starts one base BEFORE the slice, and the whole slice is adapter.
    left, right = "ACCCTGTTAT", "GAGCATGACC"   # the v11 genome track's flanks
    genome = TrackConfig("genome", "cassette", 101, 120, left, right, True)
    clip = revcomp(MOSAIC_END[1:] + I5_TAIL)
    seq = clip + right + "TTTCCAAACT"
    assert seq[len(clip)] == revcomp(MOSAIC_END[0])  # the adapter base aligned as flank
    assert extract_barcode(record(seq, f"{len(clip)}S20M", 120, 83), genome) is None


def test_genuine_reverse_read_is_kept():
    seq = LEFT + BC + RIGHT
    result = extract_barcode(record(seq, "40M", 90, 83), track())
    assert result is not None and result.barcode == BC
