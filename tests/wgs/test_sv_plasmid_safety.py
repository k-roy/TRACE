"""Portable Stage 3 functional regression tests with independent synthetic BAMs.

Requires pytest, pysam and samtools on PATH. Fixtures create all data under pytest's
temporary directory; no private data, network, audit checkout or fixture download
is required. Expected read geometry, alleles and counts are fixed before calling
the implementation. Concurrent writer cases use two bounded threads.
"""
from dataclasses import replace
from pathlib import Path
import csv
import shutil
import subprocess
import sys
import threading
from concurrent.futures import ThreadPoolExecutor

import pysam
import pytest

from trace_crispr.wgs import sv_plasmid as sv

C = 'integrated_plasmid_oligo_1'
SAMTOOLS = shutil.which('samtools')


def make_read(header, name, contig, start, seq, cigar, flag=99, mate=None, mate_start=180, tlen=160, mapq=60, qual=35):
    r = pysam.AlignedSegment(header)
    r.query_name = name
    r.query_sequence = seq
    r.flag = flag
    r.reference_id = header.get_tid(contig)
    r.reference_start = start
    r.cigarstring = cigar
    r.mapping_quality = mapq
    r.next_reference_id = header.get_tid(mate or contig)
    r.next_reference_start = mate_start
    r.template_length = tlen
    r.query_qualities = pysam.qualitystring_to_array(chr(qual+33)*len(seq))
    return r

def write_bam(path, references, records):
    header_dict = pysam.AlignmentHeader.from_references(list(references),[len(v) if isinstance(v,str) else v for v in references.values()]).to_dict()
    header_dict['PG'] = [{'ID':'fixture_align', 'CL':f'fixture_align {path.parent / "ref.fa"}'}]
    header = pysam.AlignmentHeader.from_dict(header_dict)
    reads = [make_read(header, **d) for d in records]
    with pysam.AlignmentFile(str(path),'wb',header=header) as f:
        for r in sorted(reads,key=lambda r:(r.reference_id,r.reference_start)):
            f.write(r)
    pysam.index(str(path))
    return path

def read_rows(path):
    with open(path) as f:return list(csv.DictReader(f,delimiter='\t'))


@pytest.fixture
def fixture_factory(tmp_path):
    def make(*, records=None):
        root = tmp_path / ('case' + str(len(list(tmp_path.iterdir()))))
        root.mkdir()
        refs = {'chrII': 'A' * 500, 'plasmid_demo': 'A' * 500}
        fa = root / 'ref.fa'
        fa.write_text(''.join('>' + c + '\n' + s + '\n' for c, s in refs.items()))
        pysam.faidx(str(fa))
        records = records or []
        bam = write_bam(root / 'different_from_sample.bam', refs, records)
        return dict(root=root, fa=fa, bam=bam, refs=refs, records=records)
    return make


def alignment(name='r', start=95, cigar='20M', seq='A' * 20, **kwargs):
    return dict(name=name, contig=C, start=start, cigar=cigar, seq=seq, **kwargs)

def half_b_paths(f):
    blacklist = f['root'] / 'blacklist.tsv'; blacklist.write_text('plasmid_demo\t1\t500\n')
    return dict(in_bam=str(f['bam']), out_bam=str(f['root'] / 'out.bam'), chrom='chrII', donor_start=301,
                donor_end=310, blacklist_path=str(blacklist), samtools=SAMTOOLS,
                reference_path=str(f['fa']), blacklist_reference=str(f['fa']))

@pytest.mark.parametrize('start,cigar,seq', [
    (100, '20M', 'A' * 20), (80, '20M', 'A' * 20),
    (85, '5M30D15M', 'A' * 20), (85, '5M30N15M', 'A' * 20),
])
def test_junction_requires_aligned_bases_on_both_sides(tmp_path, start, cigar, seq):
    bam = write_bam(tmp_path / 'x.bam', {C: 400}, [alignment(start=start, cigar=cigar, seq=seq)])
    assert sv.detect_plasmid_integration(str(bam), junction_distance={C: 100})['upstream_count'] == 0

@pytest.mark.parametrize('flag,mapq', [(99, 0), (99, 19), (99 | 256, 60), (99 | 2048, 60), (99 | 512, 60)])
def test_integration_mapping_quality_and_primary_filters(tmp_path, flag, mapq):
    bam = write_bam(tmp_path / 'x.bam', {C: 400}, [alignment(flag=flag, mapq=mapq)])
    assert sv.detect_plasmid_integration(str(bam), junction_distance={C: 100})['upstream_count'] == 0

@pytest.mark.parametrize('kwargs', [
    {}, {'junction_distance': {}}, {'junction_distance': {C: -1}},
    {'junction_distance': {C: 0}}, {'junction_distance': {C: 201}},
    {'junction_distance': {C: 100}, 'integrated_contigs': ['absent']},
])
def test_unknown_or_invalid_junction_never_silent_zero(tmp_path, kwargs):
    bam = write_bam(tmp_path / 'x.bam', {C: 400}, [])
    with pytest.raises(ValueError):
        sv.detect_plasmid_integration(str(bam), **kwargs)

def test_no_eligible_contig_is_unmeasured(tmp_path):
    bam = write_bam(tmp_path / 'x.bam', {'chrII': 400}, [])
    with pytest.raises(ValueError, match='eligible'):
        sv.detect_plasmid_integration(str(bam), junction_distance=100)

def test_empty_but_eligible_bam_is_valid_zero(tmp_path):
    bam = write_bam(tmp_path / 'x.bam', {C: 400}, [])
    assert sv.detect_plasmid_integration(str(bam), junction_distance={C: 100}) == {'upstream_count': 0, 'downstream_count': 0}

def test_different_donor_lengths_have_different_junctions(tmp_path):
    c2 = 'integrated_plasmid_other'
    reads = [alignment(), dict(alignment(name='r2', start=125), contig=c2)]
    bam = write_bam(tmp_path / 'x.bam', {C: 400, c2: 460}, reads)
    assert sv.detect_plasmid_integration(str(bam), junction_distance={C: 100, c2: 130}) == {'upstream_count': 2, 'downstream_count': 0}

@pytest.mark.parametrize('bad', ['chrII 1 100\n', 'chrII\t0\t100\n', 'chrII\t100\t10\n', 'chrII\tbad\t100\n'])
def test_bad_blacklist_rejected_before_processing(tmp_path, bad):
    path = tmp_path / 'blacklist.tsv'; path.write_text(bad)
    with pytest.raises(ValueError):
        sv.load_blacklist(str(path))

@pytest.mark.parametrize('flags,rname,qname', [(147, 'chrIII', 'r'), (147, 'chrII', 'other'), (99, 'chrII', 'r'), (147 | 256, 'chrII', 'r'), (147 | 8, 'chrII', 'r')])
def test_invalid_pair_cannot_supply_large_deletion(flags, rname, qname):
    h = pysam.AlignmentHeader.from_references(['chrII', 'chrIII'], [5000, 5000])
    a = make_read(h, 'r', 'chrII', 100, 'A' * 20, '10M2000D10M', flag=99, mate_start=1000)
    b = make_read(h, qname, rname, 1000, 'A' * 20, '20M', flag=flags, mate_start=100)
    assert not sv.mate_pair_large_deletion(a, b)

@pytest.mark.parametrize('flag,mapq', [(256, 60), (2048, 60), (0, 0), (512, 60)])
def test_invalid_read_cannot_supply_cigar_deletion(flag, mapq):
    h = pysam.AlignmentHeader.from_references(['chrII'], [5000])
    read = make_read(h, 'r', 'chrII', 100, 'A' * 20, '10M51D10M', flag=flag, mapq=mapq)
    assert sv.within_read_large_deletions(read) == []

def test_half_b_deletion_is_not_midpoint_coverage():
    record = sv.SamRecord('gap', 99, 'chrII', 296, '5M10D15M', '=', 350, 100)
    result = sv.classify_donor_reads([record], 'chrII', 301, 310, {})
    assert result.variant_cov == 0

@pytest.mark.parametrize('record', [
    sv.SamRecord('secondary', 99 | 256, 'chrII', 296, '20M', '=', 350, 100),
    sv.SamRecord('supplementary', 99 | 2048, 'chrII', 296, '20M', '=', 350, 100),
    sv.SamRecord('unpaired', 0, 'chrII', 296, '20M', '=', 350, 100),
    sv.SamRecord('mate_missing', 99 | 8, 'chrII', 296, '20M', '=', 350, 100),
    sv.SamRecord('mate_zero', 99, 'chrII', 296, '20M', '=', 0, 100),
])
def test_half_b_invalid_record_never_genome_count(record):
    result = sv.classify_donor_reads([record], 'chrII', 301, 310, {})
    assert result.frags_genome == result.variant_cov == 0

def test_half_b_categories_are_exclusive_in_standard_mode():
    records = [sv.SamRecord('both', 99, 'chrII', 296, '20M', '=', 350, 75), sv.SamRecord('both', 147, 'chrII', 301, '10M', '=', 301, 10)]
    result = sv.classify_donor_reads(records, 'chrII', 301, 310, {})
    assert result.unique_genome == set()
    # Standard fragments are (RG, QNAME) keys; neither record declares an RG.
    assert result.unique_plasmid == {sv.FragmentKey('', 'both')}
    assert result.included == set() and result.discarded == {sv.FragmentKey('', 'both')}
    assert result.variant_cov == 0

def test_half_b_boundary_end_before_donor_is_not_overlap():
    record = sv.SamRecord('before', 99, 'chrII', 291, '10M', '=', 350, 100)
    assert sv.classify_donor_reads([record], 'chrII', 301, 310, {}).overlapped == set()

def test_half_b_missing_reference_or_wrong_blacklist_fails(fixture_factory):
    f = fixture_factory(records=[])
    blacklist = f['root'] / 'bad.tsv'; blacklist.write_text('wrong_reference\t1\t10\n')
    out = f['root'] / 'out.bam'
    with pytest.raises(ValueError, match='reference'):
        sv.clean_donor_reads(str(f['bam']), str(out), 'chrII', 301, 310, str(blacklist), samtools=SAMTOOLS, reference_path=str(f['fa']), blacklist_reference=str(f['fa']))
    assert not out.exists() and not Path(str(out) + '.target_filtering.anno').exists()

@pytest.mark.parametrize('start,anchor,expected', [(99, 1, 1), (99, 2, 0), (98, 2, 1)])
def test_min_anchor_is_aligned_bases_per_side(tmp_path, start, anchor, expected):
    bam = write_bam(tmp_path / 'x.bam', {C: 400}, [alignment(start=start)])
    assert sv.detect_plasmid_integration(str(bam), junction_distance=100, min_anchor=anchor)['upstream_count'] == expected

def test_unknown_mapq_255_never_counts(tmp_path):
    bam = write_bam(tmp_path / 'x.bam', {C: 400}, [alignment(mapq=255)])
    assert sv.detect_plasmid_integration(str(bam), junction_distance=100)['upstream_count'] == 0

def test_valid_reciprocal_nested_pair_and_deletion_thresholds():
    h = pysam.AlignmentHeader.from_references(['chrII'], [5000])
    a = make_read(h, 'r', 'chrII', 100, 'A' * 20, '10M2000D10M', flag=99, mate_start=1000)
    b = make_read(h, 'r', 'chrII', 1000, 'A' * 20, '20M', flag=147, mate_start=100)
    assert sv.mate_pair_large_deletion(a, b)
    assert [(d.ref_start, d.ref_end, d.size) for d in sv.within_read_large_deletions(a)] == [(110, 2110, 2000)]
    a.cigarstring = '10M50D10M'
    assert sv.within_read_large_deletions(a) == []
    a.cigarstring = '10M2000N10M'
    assert not sv.mate_pair_large_deletion(a, b)
    a.cigarstring = '10M51N10M'
    assert sv.within_read_large_deletions(a) == []

@pytest.mark.parametrize('cigar', ['10M51D', '51D10M'])
def test_terminal_deletion_requires_aligned_bases_both_sides(cigar):
    h = pysam.AlignmentHeader.from_references(['chrII'], [5000])
    read = make_read(h, 'r', 'chrII', 100, 'A' * 10, cigar)
    assert sv.within_read_large_deletions(read) == []

@pytest.mark.parametrize('cigar,expected', [('5=1X14=', 1), ('5M10N15M', 0), ('5M10D15M', 0)])
def test_half_b_cigar_semantics_at_midpoint(cigar, expected):
    result = sv.classify_donor_reads([sv.SamRecord('r', 99, 'chrII', 296, cigar, '=', 350, 100)], 'chrII', 301, 310, {})
    assert result.variant_cov == expected

def test_half_b_blacklist_includes_boundary_coordinates():
    record = sv.SamRecord('r', 99, 'chrII', 296, '20M', 'plasmid', 1, 0)
    result = sv.classify_donor_reads([record], 'chrII', 301, 310, {'plasmid': [(1, 100)]})
    assert result.unique_plasmid == {sv.FragmentKey('', 'r')}

def test_half_b_parity_preserves_historical_overlapping_categories():
    records = [sv.SamRecord('both', 99, 'chrII', 296, '20M', '=', 350, 75), sv.SamRecord('both', 147, 'chrII', 301, '10M', '=', 301, 10)]
    result = sv.classify_donor_reads(records, 'chrII', 301, 310, {}, dup_mode='parity')
    assert result.unique_genome == result.unique_plasmid == {'both'}

def test_half_b_empty_indexed_bam_writes_valid_empty_outputs(fixture_factory):
    f = fixture_factory(records=[]); args = half_b_paths(f)
    result = sv.clean_donor_reads(**args)
    assert result.anno_row() == '0\t0\t0\t0'
    with pysam.AlignmentFile(args['out_bam']) as bam:
        assert list(bam) == []
    assert Path(args['out_bam'] + '.target_filtering.anno').read_text().splitlines()[1] == result.anno_row()

def test_half_b_filtered_target_records_do_not_enter_counts_or_cleaned_bam(fixture_factory):
    records = [dict(name=name, contig='chrII', start=295, seq='A' * 20, cigar='20M', mate_start=350, tlen=75, flag=flag, mapq=mapq)
               for name, flag, mapq in [('good', 99, 60), ('dup', 99 | 1024, 60), ('low', 99, 19), ('secondary', 99 | 256, 60), ('supp', 99 | 2048, 60), ('unknown', 99, 255)]]
    f = fixture_factory(records=records); args = half_b_paths(f)
    result = sv.clean_donor_reads(**args)
    assert result.anno_row() == '1\t1\t0\t0'
    with pysam.AlignmentFile(args['out_bam']) as bam:
        assert {r.query_name for r in bam} == {'good', 'dup'}

def test_half_b_conversion_failure_publishes_neither_output(fixture_factory, monkeypatch):
    f = fixture_factory(records=[]); args = half_b_paths(f)
    def failed(*a, **kw):
        raise subprocess.CalledProcessError(7, a[0])
    monkeypatch.setattr(sv.subprocess, 'run', failed)
    with pytest.raises(subprocess.CalledProcessError):
        sv.clean_donor_reads(**args)
    assert not Path(args['out_bam']).exists()
    assert not Path(args['out_bam'] + '.target_filtering.anno').exists()
    assert not list(f['root'].glob('.donor_clean_*'))

def test_half_b_same_named_contig_different_sequence_refused(fixture_factory):
    f = fixture_factory(records=[]); args = half_b_paths(f)
    wrong = f['root'] / 'different_reference.fa'
    wrong.write_text(f['fa'].read_text().replace('A', 'C')); pysam.faidx(str(wrong))
    args['blacklist_reference'] = str(wrong)
    with pytest.raises(ValueError, match='sequence identity'):
        sv.clean_donor_reads(**args)
    assert not Path(args['out_bam']).exists()

def test_half_b_renamed_identical_reference_allowed(fixture_factory):
    f = fixture_factory(records=[]); args = half_b_paths(f)
    renamed = f['root'] / 'renamed.fa'; shutil.copyfile(f['fa'], renamed); pysam.faidx(str(renamed))
    args['blacklist_reference'] = str(renamed)
    assert sv.clean_donor_reads(**args).anno_row() == '0\t0\t0\t0'

def test_half_b_source_reference_is_required(fixture_factory):
    f = fixture_factory(records=[]); args = half_b_paths(f); args.pop('blacklist_reference')
    with pytest.raises(ValueError, match='provenance'):
        sv.clean_donor_reads(**args)

@pytest.mark.parametrize('donor,extension', [('', 90), ('ACGT', 0), ('ACGT', -1)])
def test_junction_derivation_rejects_missing_sequence_or_arms(donor, extension):
    with pytest.raises(ValueError):
        sv.integrated_plasmid_junction(donor, extension)

def test_direct_summary_distinguishes_negative_and_unmeasured(tmp_path):
    positive = write_bam(tmp_path / 'positive.bam', {C: 400}, [alignment()])
    negative = write_bam(tmp_path / 'negative.bam', {C: 400}, [])
    unmeasured = write_bam(tmp_path / 'unmeasured.bam', {'chrII': 400}, [])
    rows = sv.summarize_plasmid_integration({'unmeasured': str(unmeasured), 'negative': str(negative), 'positive': str(positive)}, junction_distance=100)
    assert [(r['sample_name'], r['total_count'], r['half_a_status']) for r in rows] == [
        ('positive', 1, 'run'), ('negative', 0, 'run'), ('unmeasured', None, 'not_run_no_integrated_contigs')]

def test_half_b_cli_same_real_samtools_path(fixture_factory):
    f = fixture_factory(records=[]); args = half_b_paths(f)
    command = [sys.executable, '-m', 'trace_crispr.wgs.sv_plasmid', '--in', args['in_bam'], '--out', args['out_bam'],
               '--chr', args['chrom'], '--d-start', '301', '--d-end', '310', '--blacklist', args['blacklist_path'],
               '--reference', args['reference_path'], '--blacklist-reference', args['blacklist_reference'], '--samtools', SAMTOOLS]
    run = subprocess.run(command, capture_output=True, text=True)
    assert run.returncode == 0, run.stderr
    assert Path(args['out_bam'] + '.target_filtering.anno').read_text().splitlines()[1] == '0\t0\t0\t0'
    with pysam.AlignmentFile(args['out_bam']) as bam:
        assert list(bam) == []

def test_cleaned_bam_index_supports_random_access_and_preserves_off_target_records(fixture_factory):
    reads = [dict(name='template', contig='chrII', start=300, seq='A' * 10, cigar='10M', mate_start=300, tlen=10),
             dict(name='genome', contig='chrII', start=295, seq='A' * 20, cigar='20M', mate_start=350, tlen=75),
             dict(name='outside', contig='chrII', start=400, seq='A' * 20, cigar='20M', mate_start=450, tlen=75),
             dict(name='other_contig', contig='plasmid_demo', start=400, seq='A' * 20, cigar='20M', mate_start=450, tlen=75)]
    f = fixture_factory(records=reads); args = half_b_paths(f)
    with pysam.AlignmentFile(str(f['bam'])) as original:
        unchanged = {r.query_name: r.to_string() for r in original if r.query_name in ('outside', 'other_contig')}
    sv.clean_donor_reads(**args)
    assert Path(args['out_bam'] + '.bai').is_file()
    with pysam.AlignmentFile(args['out_bam']) as cleaned:
        assert cleaned.check_index()
        assert [r.query_name for r in cleaned.fetch('chrII', 294, 315)] == ['genome']
        kept = {r.query_name: r.to_string() for r in cleaned.fetch('chrII', 399, 425)}
        kept.update({r.query_name: r.to_string() for r in cleaned.fetch('plasmid_demo', 399, 425)})
    assert kept == unchanged

def test_half_b_large_reference_uses_csi_without_allocating_large_sequence(tmp_path):
    header = pysam.AlignmentHeader.from_references(['chrII', 'large'], [500, 700000000])
    read = make_read(header, 'far', 'large', 600000000, 'A' * 20, '20M', mate_start=600000100)
    source = tmp_path / 'in.bam'
    with pysam.AlignmentFile(str(source), 'wb', header=header) as bam:
        bam.write(read)
    pysam.index('-c', str(source))
    blacklist = tmp_path / 'blacklist'; blacklist.write_text('chrII\t1\t100\n')
    output = tmp_path / 'out.bam'
    # Explicit historical mode needs no fabricated 700-Mb reference sequence.
    sv.clean_donor_reads(str(source), str(output), 'chrII', 301, 310, str(blacklist), samtools=SAMTOOLS, dup_mode='parity')
    assert Path(str(output) + '.csi').is_file() and not Path(str(output) + '.bai').exists()
    with pysam.AlignmentFile(str(output)) as bam:
        assert [r.query_name for r in bam.fetch('large', 600000000, 600000020)] == ['far']

def test_half_b_index_failure_never_publishes_bam_or_annotation(fixture_factory, monkeypatch):
    f = fixture_factory(records=[]); args = half_b_paths(f)
    def fail(*a, **kw):
        raise RuntimeError('index construction failed')
    monkeypatch.setattr(pysam, 'index', fail)
    with pytest.raises(RuntimeError, match='index construction failed'):
        sv.clean_donor_reads(**args)
    assert not Path(args['out_bam']).exists()
    assert not Path(args['out_bam'] + '.target_filtering.anno').exists()
    assert not Path(args['out_bam'] + '.bai').exists()

def test_cleaned_bundle_bytes_identical_across_output_directories(fixture_factory):
    f = fixture_factory(records=[dict(name='genome', contig='chrII', start=295, seq='A' * 20, cigar='20M', mate_start=350, tlen=75)])
    args = half_b_paths(f)
    first = f['root'] / 'direct' / 'cleaned.bam'; second = f['root'] / 'snakemake' / 'cleaned.bam'
    args['out_bam'] = str(first); one = sv.clean_donor_reads(**args)
    args['out_bam'] = str(second); two = sv.clean_donor_reads(**args)
    assert one.anno_row() == two.anno_row() == '1\t1\t0\t0'
    for suffix in ('', '.bai', '.target_filtering.anno'):
        assert Path(str(first) + suffix).read_bytes() == Path(str(second) + suffix).read_bytes()
    with pysam.AlignmentFile(str(first)) as cleaned, pysam.AlignmentFile(str(f['bam'])) as source:
        assert cleaned.header.to_dict() == source.header.to_dict()

def test_same_output_concurrent_cleaners_never_overlap_classification(fixture_factory, monkeypatch):
    f = fixture_factory(records=[]); args = half_b_paths(f)
    original = sv.classify_donor_reads
    first_inside = threading.Event(); second_inside = threading.Event(); release = threading.Event()
    state = {'calls': 0}
    def classify(*a, **kw):
        state['calls'] += 1
        (first_inside if state['calls'] == 1 else second_inside).set()
        assert release.wait(3)
        return original(*a, **kw)
    monkeypatch.setattr(sv, 'classify_donor_reads', classify)
    with ThreadPoolExecutor(max_workers=2) as pool:
        one = pool.submit(sv.clean_donor_reads, **args)
        assert first_inside.wait(3)
        two = pool.submit(sv.clean_donor_reads, **args)
        overlapped = second_inside.wait(0.25)
        release.set(); one.result(timeout=3); two.result(timeout=3)
    assert not overlapped, 'two writers classified the same destination simultaneously'
    assert state['calls'] == 2  # standalone TRACE has serialization, not a result cache

@pytest.mark.parametrize('qname', ['', '*'])
def test_unknown_query_name_cannot_define_a_fragment(qname):
    record = sv.SamRecord(qname, 99, 'chrII', 296, '20M', '=', 350, 75)
    with pytest.raises(ValueError, match='query name'):
        sv.classify_donor_reads([record], 'chrII', 301, 310, {})

@pytest.mark.parametrize('record,expected', [
    (sv.SamRecord('left', 99, 'chrII', 301, '5S20M', '=', 350, 100), 'chimeric'),
    (sv.SamRecord('right', 99, 'chrII', 291, '20M5S', '=', 350, 100), 'chimeric'),
    (sv.SamRecord('hard_soft', 99, 'chrII', 301, '4H5S20M', '=', 350, 100), 'chimeric'),
    (sv.SamRecord('distal', 99, 'chrII', 290, '30M', 'cassette', 5, 0), 'chimeric'),
    (sv.SamRecord('template', 99, 'chrII', 301, '10M', 'cassette', 5, 0), 'plasmid'),
    (sv.SamRecord('genome', 99, 'chrII', 296, '20M', 'chrIII', 5, 0), 'genome'),
])
def test_half_b_independent_chimeric_and_distal_geometries(record, expected):
    result = sv.classify_donor_reads([record], 'chrII', 301, 310, {'cassette': [(1, 10)]})
    assert (result.frags_genome, result.frags_plasmid, result.frags_chimeric) == {
        'genome': (1, 0, 0), 'plasmid': (0, 1, 0), 'chimeric': (0, 0, 1)}[expected]

def test_fragment_length_inconsistent_with_mate_is_unresolved():
    record = sv.SamRecord('bad_tlen', 99, 'chrII', 296, '20M', '=', 350, 10)
    result = sv.classify_donor_reads([record], 'chrII', 301, 310, {})
    assert result.discarded == {sv.FragmentKey('', 'bad_tlen')}
    assert result.anno_row() == '0\t0\t0\t0'

def test_reference_builder_sequence_order_and_minus_strand_are_explicit():
    assert sv.build_plasmid_contig('AC', 'GG', 'TT', 'CA') == 'ACGGTTCA'
    assert sv.build_shifted_plasmid_contig('AACCGGTT', 4) == 'CCGGTTAA'
    assert sv.build_integrated_plasmid_contig('AAAACGTTTT', 'CG', 'AAGC', '+', bp_to_extend=4) == 'AAAACGAAGCCGTTTT'
    assert sv.build_integrated_plasmid_contig('AAAACGTTTT', 'CG', 'AAGC', '-', bp_to_extend=4) == 'AAAACGGCTTCGTTTT'
    assert sv.integrated_plasmid_junction('CG', 4) == 6

@pytest.mark.parametrize('donor,extension,strand', [('CG', 4, '?'), ('', 4, '+'), ('CG', 0, '+'), ('CG', 6, '+')])
def test_reference_builder_rejects_invalid_layout(donor, extension, strand):
    with pytest.raises(ValueError):
        sv.build_integrated_plasmid_contig('AAAACGTTTT', donor, 'AAGC', strand, bp_to_extend=extension)

def test_placeholder_query_name_cannot_link_deletion_mates():
    header = pysam.AlignmentHeader.from_references(['chrII'], [5000])
    outer = make_read(header, '*', 'chrII', 100, 'A' * 20, '10M2000D10M', flag=99, mate_start=1000)
    inner = make_read(header, '*', 'chrII', 1000, 'A' * 20, '20M', flag=147, mate_start=100)
    assert not sv.mate_pair_large_deletion(outer, inner)
