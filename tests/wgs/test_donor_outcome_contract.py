"""Portable independent donor caller contracts; expectations are authored molecules.

Run from the TRACE checkout with its test dependencies installed:
    python -m pytest tests/wgs/test_donor_outcome_contract.py
No MAGESTIC or planning-directory import is required.
"""
from pathlib import Path
import gzip
import pytest
from trace_crispr.wgs import donor_outcome as d


def strand_spec(tmp_path, *, wt='ACGT', design='CCGT', empirical=None, strand=None):
    # Coordinates and edited molecules below are authored independently of callers.
    sheet = tmp_path/'strand.tsv'
    header = 'sample\tchrom\tdonor_start_coord\tdonor_end_coord\twt_seq\tdesign_donor\tsynthetic_donor'
    row = f's\t{CHROM}\t1\t4\t{wt}\t{design}\t{empirical or design}'
    if strand is not None:
        header += '\tdonor_strand'; row += '\t'+strand
    sheet.write_text(header+'\n'+row+'\n')
    return d.read_samplesheet(sheet)[0]


@pytest.mark.parametrize('strand,pos,ref,alt', [('+',1,'A','C'), ('-',4,'T','G')])
def test_palindromic_wt_explicit_strand_selects_independently_expected_edit(tmp_path,strand,pos,ref,alt):
    spec = strand_spec(tmp_path, strand=strand)
    v,g = tmp_path/'s.vcf',tmp_path/'s.g.vcf'
    v.write_text(HEADER+vcf_line(pos,ref,alt));g.write_text(HEADER)
    result = d.precall_sample(spec,v,g,'ACGT')
    assert result.editing_outcome == 'designed_variant_installed'
    assert result.designed_variants == [f'{CHROM}_{pos}_{ref}_to_{alt}']
    assert spec.donor_strand == strand  # orientation must not mutate source metadata


@pytest.mark.parametrize('strand,designed,synthetic', [
    ('+', (1,'A','C'), (2,'C','T')), ('-', (4,'T','G'), (3,'G','A')),
])
def test_declared_strand_orients_empirical_donor_with_design(tmp_path,strand,designed,synthetic):
    spec = strand_spec(tmp_path,design='CCGT',empirical='CTGT',strand=strand)
    v,g=tmp_path/'v',tmp_path/'g'
    v.write_text(HEADER+vcf_line(*designed)+vcf_line(*synthetic));g.write_text(HEADER)
    result=d.precall_sample(spec,v,g,'ACGT')
    assert result.editing_outcome == 'designed_variant_installed;synthetic_error_installed'
    pos,ref,alt=designed
    assert result.designed_variants == [f'{CHROM}_{pos}_{ref}_to_{alt}']
    pos,ref,alt=synthetic
    assert result.synthetic_errors == [f'{CHROM}_{pos}_{ref}_to_{alt}']
    assert spec.donor_step2 == 'CTGT' and spec.donor_strand == strand


@pytest.mark.parametrize('design,empirical', [('CCGT','CCGT'),('ACGT','CCGT')])
def test_palindromic_wt_requires_metadata_for_either_asymmetric_donor(tmp_path,design,empirical):
    spec = strand_spec(tmp_path,design=design,empirical=empirical)
    with pytest.raises(ValueError,match='donor_strand'):
        d.orient_to_genome(spec,'ACGT')


@pytest.mark.parametrize('strand',['plus','1','.'])
def test_invalid_donor_strand_refused(tmp_path,strand):
    with pytest.raises(ValueError,match='donor_strand'):
        strand_spec(tmp_path,strand=strand)


@pytest.mark.parametrize('wt,strand',[('ACGA','-'),('TCGT','+')])
def test_declared_donor_strand_must_match_reference(tmp_path,wt,strand):
    spec = strand_spec(tmp_path,wt=wt,design='CCGA',strand=strand)
    with pytest.raises(ValueError,match='donor_strand'):
        d.orient_to_genome(spec,'ACGA')


def test_symmetric_wt_design_and_empirical_need_no_strand_metadata(tmp_path):
    spec = strand_spec(tmp_path,design='ACGT',empirical='ACGT')
    oriented = d.orient_to_genome(spec,'ACGT')
    assert oriented.donor_wt == oriented.donor_design == oriented.donor_step2 == 'ACGT'


@pytest.mark.parametrize('source', ['vcf', 'gvcf'])
def test_symbolic_event_starting_before_donor_retains_overlap_evidence(tmp_path, source):
    event = vcf_line(10, ref='N', alt='<DEL>', fmt='GT', value='1', info='END=15')
    block = vcf_line(11, alt='<NON_REF>', fmt='GT:MIN_DP', value='0:8', info='END=20')
    kwargs = {'vcf': event, 'gvcf': block} if source == 'vcf' else {'gvcf': event}
    result = call(tmp_path, **kwargs)
    assert result.editing_outcome == 'unknown;unresolved_allele'
    assert result.ambiguous_genotypes == [
        {'chrom': CHROM, 'pos': 10, 'GT': '1', 'ALT': '<DEL>', 'reason': 'unresolved_allele'}
    ]


def test_minus_strand_donor_and_gzipped_variant_inputs(tmp_path):
    v, g = tmp_path/'v.vcf.gz', tmp_path/'g.vcf.gz'
    with gzip.open(v, 'wt') as fh: fh.write(HEADER+vcf_line(20))
    with gzip.open(g, 'wt') as fh: fh.write(HEADER)
    spec = d.DonorSpec('s', CHROM, 11, 20, d.reverse_complement(WT),
                       d.reverse_complement(WT[:-1]+'C'), d.reverse_complement(WT[:-1]+'C'))
    result = d.precall_sample(spec, v, g, REF)
    assert result.editing_outcome == 'designed_variant_installed'
    assert result.designed_variants == [f'{CHROM}_20_A_to_C']


@pytest.mark.parametrize('info', ['END', 'END=abc'])
def test_invalid_end_value_is_explicit_input_error(tmp_path, info):
    with pytest.raises(ValueError, match='END'):
        call(tmp_path, gvcf=vcf_line(11, alt='<NON_REF>', fmt='GT:MIN_DP', value='0:8', info=info))

HEADER = '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts\n'
CHROM = 'integrated_plasmid_oligo_1'
WT = 'ACGTACGTAA'
REF = 'N' * 10 + WT + 'N' * 10


def vcf_line(pos, ref='A', alt='C', fmt='GT:AD', value='1:0,8', info='.', chrom=CHROM):
    return f'{chrom}\t{pos}\t.\t{ref}\t{alt}\t50\tPASS\t{info}\t{fmt}\t{value}\n'


def call(tmp_path, vcf='', gvcf='', design='ACGTACGTAC', **kwargs):
    vpath, gpath = tmp_path/'s.vcf', tmp_path/'s.g.vcf'
    vpath.write_text(HEADER + vcf)
    gpath.write_text(HEADER + gvcf)
    spec = d.DonorSpec('s', CHROM, 11, 20, WT, design, design)
    return d.precall_sample(spec, vpath, gpath, REF, **kwargs)


@pytest.mark.parametrize('position', range(11, 21))
def test_each_donor_position_including_terminal_and_underscore(tmp_path, position):
    i = position - 11
    alt = 'C' if WT[i] != 'C' else 'A'
    design = WT[:i] + alt + WT[i+1:]
    result = call(tmp_path, vcf_line(position, WT[i], alt), design=design)
    assert result.editing_outcome == 'designed_variant_installed'
    assert result.designed_variants == [f'{CHROM}_{position}_{WT[i]}_to_{alt}']


@pytest.mark.parametrize('ref,alt,design', [('A','AC', WT+'C'), ('AA','A',WT[:-1])])
def test_terminal_insertion_or_deletion(tmp_path, ref, alt, design):
    pos = 21 - len(ref)
    assert call(tmp_path, vcf_line(pos,ref,alt), design=design).editing_outcome == 'designed_variant_installed'


def test_two_nonoverlapping_indels_reconstruct_in_reference_coordinates(tmp_path):
    # Independently write the intended molecule: insert T after local A0, delete A8.
    design = 'ATCGTACGTA'
    records = vcf_line(11,'A','AT') + vcf_line(18,'TA','T')
    assert call(tmp_path, records, design=design).editing_outcome == 'designed_variant_installed'


def test_ref_mismatch_fails_instead_of_mutating_other_reference(tmp_path):
    with pytest.raises(ValueError, match='REF'):
        call(tmp_path, vcf_line(20,'G','C'))


@pytest.mark.parametrize('fmt,value', [('GT:AD','1:0,8'), ('AD:DP:GT','0,8:8:1'), ('GT','1')])
def test_genotype_and_format_order(tmp_path, fmt, value):
    assert call(tmp_path, vcf_line(20,fmt=fmt,value=value)).editing_outcome == 'designed_variant_installed'


def test_multiallelic_selects_genotype_not_first_alt(tmp_path):
    result = call(tmp_path, vcf_line(20,alt='G,C',fmt='AD:GT',value='0,0,8:2'))
    assert result.editing_outcome == 'designed_variant_installed'
    assert result.designed_variants == [f'{CHROM}_20_A_to_C']


@pytest.mark.parametrize('gt', ['3','-1','x'])
def test_ambiguous_or_invalid_genotype_refused(tmp_path, gt):
    with pytest.raises(ValueError):
        call(tmp_path, vcf_line(20,alt='G,C',fmt='GT',value=gt))


def test_reference_genotype_does_not_inherit_plate_variant(tmp_path):
    assert call(tmp_path, vcf_line(20,value='0:8,0')).editing_outcome == 'unknown'


def test_missing_genotype_does_not_invent_call(tmp_path):
    assert call(tmp_path, vcf_line(20,value='./.:0,8')).editing_outcome == 'unknown;ambiguous_genotype'


def test_exact_duplicate_deduplicated_conflict_refused(tmp_path):
    line = vcf_line(20)
    assert len(call(tmp_path,line+line).designed_variants) == 1
    with pytest.raises(ValueError, match='overlapping'):
        call(tmp_path,line+vcf_line(20,alt='G'))


def test_low_confidence_terminal_edit_uses_structured_id(tmp_path):
    result = call(tmp_path, gvcf=vcf_line(20))
    assert result.editing_outcome == 'designed_variant_installed'
    assert result.designed_variants == [f'(LC){CHROM}_20_A_to_C']


def test_homref_gvcf_alt_list_does_not_promote_uncalled_alt(tmp_path):
    assert call(tmp_path,gvcf=vcf_line(20,value='0:8,0')).editing_outcome == 'unknown'


def test_min_dp_end_block_covers_whole_window(tmp_path):
    block=vcf_line(11,alt='<NON_REF>',fmt='MIN_DP:GT:DP',value='8:0:10',info='END=20')
    assert call(tmp_path,gvcf=block).editing_outcome == 'designed_variant_unedited'


def test_ad_only_end_block_does_not_claim_minimum_interior_depth(tmp_path):
    block=vcf_line(11,alt='<NON_REF>',value='0:8,0',info='END=20')
    assert call(tmp_path,gvcf=block).editing_outcome == 'unknown'


def test_repeated_single_position_cannot_fill_window(tmp_path):
    repeated=vcf_line(11,alt='<NON_REF>',value='0:8,0') * 10
    assert call(tmp_path,gvcf=repeated).editing_outcome == 'unknown'


@pytest.mark.parametrize('fmt,value', [('GT:AD','0:8,0'), ('GT:PL:DP:AD','0:0,9:8:8,0'),('DP:GT','8:0')])
def test_explicit_reference_positions_use_named_depth(tmp_path, fmt, value):
    lines=''.join(vcf_line(p,WT[p-11],'<NON_REF>',fmt,value) for p in range(11,21))
    assert call(tmp_path,gvcf=lines).editing_outcome == 'designed_variant_unedited'


@pytest.mark.parametrize('kwargs', [{'ne_threshold':0},{'quantification_window_size':0}, {'quantification_window_cov_perc':0},{'quantification_window_cov_perc':1.1}])
def test_invalid_coverage_parameters_fail(tmp_path, kwargs):
    with pytest.raises(ValueError): call(tmp_path, **kwargs)


def test_short_donor_uses_its_actual_window_denominator(tmp_path):
    spec=d.DonorSpec('s','chr1',1,2,'AA','AC','AC')
    v,g=tmp_path/'v',tmp_path/'g';v.write_text(HEADER)
    g.write_text(HEADER+vcf_line(1,alt='<NON_REF>',fmt='GT:MIN_DP',value='0:5',info='END=2',chrom='chr1'))
    assert d.precall_sample(spec,v,g,'AA').editing_outcome == 'designed_variant_unedited'


def test_loaded_reference_rejects_different_or_changed_source(tmp_path):
    original = tmp_path/'original.fa'; original.write_text('>chr1\nAAAA\n')
    loaded = d.load_reference(original)
    assert loaded.require_current(original) == {'chr1': 'AAAA'}
    other = tmp_path/'other.fa'; other.write_text(original.read_text())
    with pytest.raises(ValueError, match='reference'):
        loaded.require_current(other)
    original.write_text('>chr1\nAAAAA\n')
    with pytest.raises(ValueError, match='reference'):
        loaded.require_current(original)


def test_reference_identical_design_is_not_called_installed_without_data(tmp_path):
    assert call(tmp_path,design=WT).editing_outcome == 'unknown'


@pytest.mark.parametrize('gvcf', [
    vcf_line(20,'T','<NON_REF>',value='0:8,0'),
    vcf_line(11,alt='<NON_REF>',fmt='GT:MIN_DP',value='0:8',info='END=100'),
    vcf_line(11,alt='<NON_REF>',fmt='GT:MIN_DP',value='0:8',info='END=10'),
    vcf_line(11,alt='C',fmt='GT:MIN_DP',value='1:8',info='END=20'),
    vcf_line(11,alt='<NON_REF>',value='0:-1,0'),
])
def test_invalid_reference_coverage_evidence_refused(tmp_path,gvcf):
    with pytest.raises(ValueError):call(tmp_path,gvcf=gvcf)


@pytest.mark.parametrize('gt', ['0/1','1/.','1/2','0|1','./.'])
def test_ambiguous_genotype_is_explicit_unknown_with_raw_evidence(tmp_path,gt):
    result=call(tmp_path,vcf_line(20,alt='C,G',fmt='GT',value=gt))
    assert result.editing_outcome=='unknown;ambiguous_genotype'
    assert result.ambiguous_genotypes==[{'chrom':CHROM,'pos':20,'GT':gt}]
    assert result.designed_variants==[]


@pytest.mark.parametrize('gt', ['1','1/1','1|1'])
def test_single_supported_genotype_is_callable(tmp_path,gt):
    assert call(tmp_path,vcf_line(20,fmt='GT',value=gt)).editing_outcome=='designed_variant_installed'


def test_homozygous_reference_needs_coverage_then_unedited(tmp_path):
    vcf=vcf_line(20,fmt='GT',value='0/0')
    assert call(tmp_path,vcf).editing_outcome=='unknown'
    block=vcf_line(11,alt='<NON_REF>',fmt='GT:MIN_DP',value='0/0:8',info='END=20')
    assert call(tmp_path,vcf,gvcf=block).editing_outcome=='designed_variant_unedited'


def test_donor_crossing_variant_is_not_silently_ignored(tmp_path):
    with pytest.raises(ValueError,match='boundary'):
        call(tmp_path,vcf_line(19,'AAN','A'))


@pytest.mark.parametrize('content', ['', vcf_line(20), HEADER+HEADER, HEADER+vcf_line(20).rstrip()+'\textra\n'])
def test_vcf_header_and_width_contract(tmp_path,content):
    v,g=tmp_path/'v',tmp_path/'g';v.write_text(content);g.write_text(HEADER)
    spec=d.DonorSpec('s',CHROM,11,20,WT,WT[:-1]+'C',WT[:-1]+'C')
    with pytest.raises(ValueError,match='header|column'):
        d.precall_sample(spec,v,g,REF)


@pytest.mark.parametrize('text', ['ACGT\n>chr\nACGT\n', '>chr\n\n', '>chr\nACGT\n>chr\nTGCA\n', '>chr\nAC*GT\n'])
def test_malformed_reference_fasta_is_not_silently_normalized(tmp_path,text):
    path=tmp_path/'bad.fa';path.write_text(text)
    with pytest.raises(ValueError):d.load_fasta(str(path))


def long_donor_call(tmp_path, design, gvcf):
    v,g=tmp_path/'long.vcf',tmp_path/'long.gvcf'
    v.write_text(HEADER);g.write_text(HEADER+gvcf)
    return d.precall_sample(d.DonorSpec('s',CHROM,1,20,'A'*20,design,design),v,g,'A'*20)


def test_unedited_requires_observing_intended_position_outside_central_window(tmp_path):
    center=vcf_line(6,alt='<NON_REF>',fmt='GT:MIN_DP',value='0:8',info='END=15')
    assert long_donor_call(tmp_path,'C'+'A'*19,center).editing_outcome=='unknown'


def test_unedited_observed_intended_position_plus_central_window(tmp_path):
    center=vcf_line(6,alt='<NON_REF>',fmt='GT:MIN_DP',value='0:8',info='END=15')
    diagnostic=vcf_line(1,alt='<NON_REF>',value='0:8,0')
    assert long_donor_call(tmp_path,'C'+'A'*19,center+diagnostic).editing_outcome=='designed_variant_unedited'


@pytest.mark.parametrize('design',['A'*21,'A'*19])
def test_unedited_indel_design_requires_complete_donor_reference_evidence(tmp_path,design):
    center=vcf_line(6,alt='<NON_REF>',fmt='GT:MIN_DP',value='0:8',info='END=15')
    assert long_donor_call(tmp_path,design,center).editing_outcome=='unknown'
    complete=vcf_line(1,alt='<NON_REF>',fmt='GT:MIN_DP',value='0:8',info='END=20')
    assert long_donor_call(tmp_path,design,complete).editing_outcome=='designed_variant_unedited'


@pytest.mark.parametrize('alt',['<NON_REF>','<*>','<DEL>','*'])
def test_selected_symbolic_allele_is_explicit_unknown(tmp_path,alt):
    records=vcf_line(11,alt='<NON_REF>',fmt='GT:MIN_DP',value='0:8',info='END=19')
    records+=vcf_line(20,alt=alt,fmt='GT:AD',value='1:8,8')
    outcome=call(tmp_path,gvcf=records)
    assert outcome.editing_outcome=='unknown;unresolved_allele'
    assert outcome.ambiguous_genotypes==[{'chrom':CHROM,'pos':20,'GT':'1','ALT':alt,'reason':'unresolved_allele'}]


def test_gt_cannot_select_missing_alt(tmp_path):
    with pytest.raises(ValueError,match='missing ALT'):
        call(tmp_path,vcf_line(20,alt='.',fmt='GT',value='1'))


@pytest.mark.parametrize('record', [
    vcf_line(11,alt='<NON_REF>',fmt='GT:MIN_DP',value='0:8',info='END=15;END=20'),
    vcf_line(20,fmt='GT::AD',value='1:.:0,8'),
    vcf_line(20,alt='C,G',value='2:8,0'),
    vcf_line(20,alt='C',value='1:8,-1'),
    vcf_line(20,alt='C,C',value='1:0,8,0'),
])
def test_malformed_gvcf_field_semantics_refused(tmp_path,record):
    with pytest.raises(ValueError):call(tmp_path,gvcf=record)


def test_trace_requires_literal_format_header(tmp_path):
    v,g=tmp_path/'bad_header.vcf',tmp_path/'g'
    v.write_text(HEADER.replace('\tFORMAT\t','\tOTHER\t')+vcf_line(20));g.write_text(HEADER)
    spec=d.DonorSpec('s',CHROM,11,20,WT,WT[:-1]+'C',WT[:-1]+'C')
    with pytest.raises(ValueError,match='header'):d.precall_sample(spec,v,g,REF)


@pytest.mark.parametrize('alt',['<DEL>','<NON_REF>'])
def test_selected_symbolic_end_is_not_a_reference_block(tmp_path,alt):
    record=vcf_line(15,alt=alt,fmt='GT:MIN_DP',value='1:8',info='END=20')
    outcome=call(tmp_path,gvcf=record)
    assert outcome.editing_outcome=='unknown;unresolved_allele'


def test_conflicting_gvcf_blocks_refused_but_exact_duplicates_idempotent(tmp_path):
    first=vcf_line(11,alt='<NON_REF>',fmt='GT:MIN_DP',value='0:8',info='END=20')
    second=vcf_line(11,alt='<NON_REF>',fmt='GT:MIN_DP',value='0:0',info='END=20')
    with pytest.raises(ValueError,match='overlapping'):
        call(tmp_path,gvcf=first+second)
    assert call(tmp_path,gvcf=first+first).editing_outcome=='designed_variant_unedited'
