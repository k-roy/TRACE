"""A no-call keeps the read counters that explain it (audit W06).

At 9d654b8 the no-call emit wrote a blank Bc1Call, erasing n_reads_seen and n_contaminant,
so "no reads here" and "every read here was dropped as contamination" looked identical. On
plate 5031 the kept counters are what showed s25 as 3 reads seen, 2 of them contaminant.
"""
import pysam

from trace_crispr.wgs.bc1_consensus import TrackConfig, call_sample, revcomp

BC = "ACGTCAGTTCGATGACCTGA"
LEFT = "GGTCATGCTC"
RIGHT = "ATAACAGGGT"


def record(name, seq, cigar, start, flag):
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


def test_no_call_keeps_the_contaminant_count(tmp_path):
    # Two barcode-bearing reads that are not proper pairs (flag 97): both are dropped.
    reads = [record(f"r{i}", LEFT + BC, "10M20S", 90, 97) for i in range(2)]
    track = TrackConfig("genome", "cassette", 101, 120, LEFT, RIGHT)
    with pysam.AlignmentFile(str(bam_file(tmp_path, reads))) as bam:
        call = call_sample(bam, "well_A1", track)
    assert (call.track, call.barcode) == ("none", "")
    assert (call.n_reads_seen, call.n_contaminant, call.n_amplicon) == (2, 2, 0)


def test_no_call_keeps_the_amplicon_count(tmp_path):
    # A proper pair with a randomer clipped at a left-window start: the amplicon filter drops it.
    reads = [record("amp", "ACGTA" + LEFT + revcomp(BC), "5S10M20S", 90, 99)]
    track = TrackConfig("genome", "cassette", 101, 120, LEFT, RIGHT, amplicon_left_window=(91, 91))
    with pysam.AlignmentFile(str(bam_file(tmp_path, reads))) as bam:
        call = call_sample(bam, "well_A1", track)
    assert (call.track, call.barcode) == ("none", "")
    assert (call.n_reads_seen, call.n_contaminant, call.n_amplicon) == (1, 1, 1)
