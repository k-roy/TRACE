"""Portable reconciliation molecular fixtures and exact expected outcomes.

All BAMs/FASTA inputs are generated here from independently specified molecules.
Requires pytest and pysam; no external datasets or audit-workspace imports.
"""
import csv
import json
from pathlib import Path
from types import SimpleNamespace

import pysam
import pytest

from trace_crispr.wgs import reconcile as r

BC = 'ACGTCAGTTCGATGACCTGA'
GLEFT, GRIGHT = 'ACCCTGTTAT', 'GAGCATGACC'
PLEFT, PRIGHT = 'GGTCATGCTC', 'ATAACAGGGT'

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

@pytest.fixture
def fixture_factory(tmp_path):
    def make(*, records=None, plasmid=True, bridge=True, window=(301,310), extra_refs=None):
        root=tmp_path/('case'+str(len(list(tmp_path.iterdir()))));root.mkdir()
        refs={'chrII':'A'*90+GLEFT+'N'*20+GRIGHT+'A'*370}
        if plasmid:refs['plasmid_demo']='A'*90+PLEFT+'N'*20+PRIGHT+'A'*370
        refs.update(extra_refs or {})
        fa=root/'ref.fa';fa.write_text(''.join('>'+c+'\n'+s+'\n' for c,s in refs.items()));pysam.faidx(str(fa))
        if records is None:
            from trace_crispr.wgs.bc1_consensus import revcomp
            records=[dict(name='g1',contig='chrII',start=90,seq=GLEFT+revcomp(BC),cigar='10M20S'),dict(name='g2',contig='chrII',start=89,seq='A'+GLEFT+revcomp(BC),cigar='11M20S',tlen=161)]
            if plasmid:records += [dict(name='p1',contig='plasmid_demo',start=90,seq=PLEFT+BC,cigar='10M20S'),dict(name='p2',contig='plasmid_demo',start=89,seq='A'+PLEFT+BC,cigar='11M20S',tlen=161)]
        bam=write_bam(root/'different_from_sample.bam',refs,records)
        manifest=root/'manifest.tsv';manifest.write_text('sample\tbam\treference\nwell_A1\t'+str(bam)+'\t'+str(fa)+'\n')
        intended='[{"pos":301,"ref":"A","alt":"C"}]' if window[0] <= 301 <= window[1] else '[]'
        windows=root/'windows.tsv';windows.write_text('oligo_name\tchrom\tdonor_start_coord\tdonor_end_coord\tintended_variants\noligo_1\tchrII\t%d\t%d\t%s\n'%(*window,intended))
        bt=root/'bridge.tsv';bt.write_text('bc1\toligo_name\n'+BC+'\toligo_1\n')
        sidecar={'barcode_locus':'chrII:101-120','barcode_left_sequence':GLEFT,'barcode_right_sequence':GRIGHT,'plasmid_bc_locus':'plasmid_demo:101-120','plasmid_bc_left_sequence':PLEFT,'plasmid_bc_right_sequence':PRIGHT}
        flat=root/'locus.txt';flat.write_text(''.join(k+": '"+v+"'\n" for k,v in sidecar.items()))
        return dict(root=root,fa=fa,bam=bam,flat=flat,refs=refs,bridge=bt,windows=windows,records=records)
    return make

def read_rows(path):
    with open(path) as f:return list(csv.DictReader(f,delimiter='\t'))


def read(name, seq, *, start=300, cigar=None, flag=0, qual=35):
    return dict(name=name, contig='chrII', start=start, seq=seq,
                cigar=cigar or f'{len(seq)}M', flag=flag, qual=qual, mapq=60)


def site(*alleles):
    return r.EditSite('chrII', 301, 320, 'oligo_1', tuple(r.IntendedAllele(*a) for a in alleles))


def score(f, es):
    with pysam.AlignmentFile(str(f['bam'])) as bam, pysam.FastaFile(str(f['fa'])) as fa:
        return r.score_edit_site(bam, fa, es)


def windows(f, alleles, *, name='oligo_1'):
    f['windows'].write_text('oligo_name\tchrom\tdonor_start_coord\tdonor_end_coord\tintended_variants\n'
        + name+'\tchrII\t301\t320\t'+json.dumps([dict(zip(('pos','ref','alt'), a)) for a in alleles])+'\n')


@pytest.mark.parametrize('alts,total,status', [(0,4,'wt'),(1,4,'partial'),(2,4,'partial'),(4,4,'full'),(3,3,'nocov'),(0,0,'nocov')])
def test_exact_snp_denominators(fixture_factory, alts, total, status):
    f=fixture_factory(records=[read(f'r{i}', ('C' if i<alts else 'A')+'A'*19) for i in range(total)])
    call=score(f,site((301,'A','C')))
    assert (call.status,call.cov,call.vaf)==(status,total,alts/total if total else 0)
    e=call.allele_evidence[0]
    assert (e.alt_count,e.ref_count,e.other_count,e.depth)==(alts,total-alts,0,total)


@pytest.mark.parametrize('seq,status,alt,ref,other', [('CAAAAAAAAAAAAAAAAAAA','full',4,0,0),
    ('ACAAAAAAAAAAAAAAAAAA','wt',0,4,0),('GAAAAAAAAAAAAAAAAAAA','other',0,0,4)])
def test_nearby_noise_and_wrong_alt_are_not_designed_edit(fixture_factory,seq,status,alt,ref,other):
    f=fixture_factory(records=[read(f'r{i}',seq) for i in range(4)])
    call=score(f,site((301,'A','C')))
    assert call.status==status
    e=call.allele_evidence[0]
    assert (e.alt_count,e.ref_count,e.other_count)==(alt,ref,other)
    cell=r.reconcile_cell('s',[call],bc_confidence='high')
    assert cell.label==('SINGLE_EDIT' if alt else 'NO_EDIT')


@pytest.mark.parametrize('allele,seq,cigar,expected', [
    ((301,'AA','CG'),'CG'+'A'*18,'20M','full'),
    ((301,'AA','CG'),'CA'+'A'*18,'20M','other'),
    ((305,'A','AC'),'A'*5+'C'+'A'*15,'5M1I15M','full'),
    ((305,'A','AC'),'A'*5+'G'+'A'*15,'5M1I15M','other'),
    ((305,'AA','A'),'A'*19,'5M1D14M','full'),
    ((305,'AAA','A'),'A'*19,'5M1D14M','other'),
    ((305,'AA','A'),'A'*19,'5M1N14M','nocov'),
])
def test_exact_mnv_and_indel_identity(fixture_factory,allele,seq,cigar,expected):
    f=fixture_factory(records=[read(f'r{i}',seq,cigar=cigar) for i in range(4)])
    call=score(f,site(allele))
    assert call.status==expected
    assert call.cov==(0 if expected=='nocov' else 4)


@pytest.mark.parametrize('flag',[256,512,1024,2048])
def test_nonprimary_duplicate_qcfail_do_not_supply_alt(fixture_factory,flag):
    f=fixture_factory(records=[read('bad','C'+'A'*19,flag=flag)]+[read(f'wt{i}','A'*20) for i in range(4)])
    e=score(f,site((301,'A','C'))).allele_evidence[0]
    assert (e.alt_count,e.depth)==(0,4)


@pytest.mark.parametrize('bad_seq,qual',[('C'+'A'*19,0),('N'+'A'*19,35)])
def test_low_quality_and_ambiguous_bases_are_missing(fixture_factory,bad_seq,qual):
    f=fixture_factory(records=[read('bad',bad_seq,qual=qual)]+[read(f'wt{i}','A'*20) for i in range(4)])
    assert score(f,site((301,'A','C'))).cov==4


@pytest.mark.parametrize('second,depth,alts',[('C',1,1),('A',0,0)])
def test_overlapping_mates_one_fragment_or_discordant_missing(fixture_factory,second,depth,alts):
    f=fixture_factory(records=[read('pair','C'+'A'*19,flag=99),read('pair',second+'A'*19,flag=147)])
    e=score(f,site((301,'A','C'))).allele_evidence[0]
    assert (e.depth,e.alt_count)==(depth,alts)


def test_multiallele_design_partial_installation_is_not_mixture(fixture_factory):
    f=fixture_factory(records=[read(f'r{i}','C'+'A'*19) for i in range(4)])
    call=score(f,site((301,'A','C'),(310,'A','G')))
    cell=r.reconcile_cell('s',[call],bc_confidence='high')
    assert (call.status,cell.label,cell.n_incomplete)==('incomplete','PARTIAL_EDIT',1)
    assert cell.confidence=='low'
    assert 'PARTIAL_VAF' not in cell.flags


def test_two_intended_alleles_in_one_oligo_count_as_one_edit(fixture_factory):
    f=fixture_factory(records=[read(f'r{i}','C'+'A'*8+'G'+'A'*10) for i in range(4)])
    call=score(f,site((301,'A','C'),(310,'A','G')))
    cell=r.reconcile_cell('s',[call],bc_confidence='high')
    assert (call.status,cell.n_full,cell.label)==('full',1,'SINGLE_EDIT')


def test_unassessed_legacy_window_never_claims_designed_edit(fixture_factory):
    f=fixture_factory(records=[read(f'r{i}','C'+'A'*19) for i in range(4)])
    call=score(f,site())
    cell=r.reconcile_cell('s',[call],bc_confidence='high')
    assert (call.vaf,call.cov,call.status)==(1,4,'unassessed')
    assert (cell.label,cell.n_full,cell.n_unassessed,cell.confidence)==('LOW_CONF',0,1,'low')


def test_covered_maximum_beats_singleton_in_legacy_window(fixture_factory):
    f=fixture_factory(records=[read('singleton','C')]+[read(f'r{i}','C' if i<9 else 'A',start=309) for i in range(10)])
    call=score(f,site())
    assert (call.vaf,call.cov,call.status)==(.9,10,'unassessed')


def test_uncovered_other_design_lowers_single_edit_confidence():
    c=r.reconcile_cell('s',[r.SiteCall(site(),1,4,'full'),r.SiteCall(site(),0,0,'nocov')],bc_confidence='high')
    assert c.label=='SINGLE_EDIT' and c.confidence=='low' and 'UNCOVERED_SITES' in c.flags


@pytest.mark.parametrize('n_bc,partial,label',[(1,False,'SINGLE_EDIT'),(2,False,'CONTAMINATION'),(1,True,'CONTAMINATION')])
def test_contamination_is_explicit_mixture_signal(n_bc,partial,label):
    c=r.reconcile_cell('s',[r.SiteCall(site(),.5 if partial else 1,4,'partial' if partial else 'full')],n_genome_bc=n_bc,bc_confidence='high')
    assert c.label==label
    assert ('MIXTURE_SIGNAL_NOT_PROOF' in c.flags)==(label=='CONTAMINATION')


def test_unknown_bc1_lowers_confidence_and_does_not_establish_absence():
    c=r.reconcile_cell('s',[r.SiteCall(site(),0,4,'wt')],bc_confidence='high',flags=('UNRESOLVED_BC1',))
    assert c.label=='NO_EDIT' and c.confidence=='low'


@pytest.mark.parametrize('kwargs', [{'vaf':float('nan'),'cov':4},{'vaf':2,'cov':4},{'vaf':0,'cov':-1},
    {'vaf':0,'cov':4,'tau_wt':.9,'tau_full':.8},{'vaf':0,'cov':4,'min_cov':0}])
def test_invalid_classification_values_fail(kwargs):
    with pytest.raises(ValueError):r.classify_site(**kwargs)


@pytest.mark.parametrize('row', ['a\tchrII\t0\t10', 'a\tchrII\t20\t10','a\tchrII\tx\t10','a\tchrII\t1',
    '\tchrII\t1\t10','a\t\t1\t10'])
def test_invalid_window_rows_fail(tmp_path,row):
    p=tmp_path/'w.tsv';p.write_text('oligo_name\tchrom\tdonor_start_coord\tdonor_end_coord\n'+row+'\n')
    with pytest.raises(ValueError):r.load_oligo_windows(str(p))


@pytest.mark.parametrize('variants',[ [{'pos':301,'ref':'A','alt':'A'}], [{'pos':300,'ref':'A','alt':'C'}],
    [{'pos':301,'ref':'A','alt':'<DEL>'}], [{'pos':301,'ref':'AA','alt':'G'}],
    [{'pos':301,'ref':'A','alt':'C'},{'pos':301,'ref':'A','alt':'G'}],
    [{'pos':'301','ref':'A','alt':'C'}], {'pos':301,'ref':'A','alt':'C'}])
def test_invalid_intended_alleles_fail(fixture_factory,variants):
    f=fixture_factory()
    f['windows'].write_text('oligo_name\tchrom\tdonor_start_coord\tdonor_end_coord\tintended_variants\na\tchrII\t301\t320\t'+json.dumps(variants)+'\n')
    with pytest.raises(ValueError):r.load_oligo_windows(str(f['windows']))


def test_reference_mismatch_fails_before_edit_claim(fixture_factory):
    f=fixture_factory(records=[read('r','C'+'A'*19)])
    with pytest.raises(ValueError,match='does not match reference'):score(f,site((301,'G','C')))


def test_conflicting_sites_and_missing_sites_rejected():
    with pytest.raises(ValueError,match='no donor window'):r.build_edit_sites(['missing'],{})
    with pytest.raises(ValueError,match='(conflicting|overlapping) intended designs'):
        r.build_edit_sites(['a','b'],{'a':site((301,'A','C')),'b':site((301,'A','G'))})


@pytest.mark.parametrize('allele,seq,cigar', [((305,'AA','A'),'A'*19,'8M1D11M'),
    ((305,'A','AA'),'A'*21,'8M1I12M')])
def test_repeat_equivalent_indel_representation_is_not_a_false_wt(fixture_factory,allele,seq,cigar):
    # Same molecular haplotype in the A repeat, but a different legal CIGAR placement.
    f=fixture_factory(records=[read(f'r{i}',seq,cigar=cigar) for i in range(4)])
    call=score(f,site(allele))
    assert (call.status,call.allele_evidence[0].alt_count,call.cov)==('full',4,4)


def test_same_intended_event_in_two_shifted_windows_is_one_edit():
    a=site((301,'A','C'))
    b=r.EditSite('chrII',300,319,'other',(r.IntendedAllele(301,'A','C'),))
    assert len(r.build_edit_sites(['a','b'],{'a':a,'b':b}))==1


def test_nested_designs_do_not_double_count_shared_intended_event():
    a=site((301,'A','C'))
    b=r.EditSite('chrII',300,319,'other',(r.IntendedAllele(301,'A','C'),r.IntendedAllele(310,'A','G')))
    with pytest.raises(ValueError,match='overlapping intended'):
        r.build_edit_sites(['a','b'],{'a':a,'b':b})


def test_reconcile_cli_zero_glob_refuses_output(fixture_factory):
    f=fixture_factory();windows(f,[(301,'A','C')]);out=f['root']/'out.tsv'
    args=r._build_arg_parser().parse_args(['--bam-dir',str(f['root']),'--bam-suffix','.absent',
        '--fasta',str(f['fa']),'--config',str(f['flat']),'--bridge-table',str(f['bridge']),
        '--oligo-window-table',str(f['windows']),'--out',str(out)])
    with pytest.raises(ValueError,match='no BAMs match'):r.run(args)
    assert not out.exists()


def test_reconcile_cli_failure_does_not_leave_success_shaped_prefix(fixture_factory,monkeypatch):
    f=fixture_factory();windows(f,[(301,'A','C')]);out=f['root']/'out.tsv'
    out.write_text('OLD COMPLETE OUTPUT\n')
    args=r._build_arg_parser().parse_args(['--bam-dir',str(f['root']),
        '--fasta',str(f['fa']),'--config',str(f['flat']),'--bridge-table',str(f['bridge']),
        '--oligo-window-table',str(f['windows']),'--out',str(out)])
    def fail(*a,**kw):raise ValueError('deliberate scoring failure')
    monkeypatch.setattr(r,'score_edit_site',fail)
    with pytest.raises(ValueError,match='deliberate'):r.run(args)
    assert out.read_text()=='OLD COMPLETE OUTPUT\n'
    assert not list(f['root'].glob('.reconcile-*'))


def test_reconcile_cli_default_suffix_matches_bam_and_counts_actual_allele(fixture_factory):
    f=fixture_factory();windows(f,[(301,'A','C')]);out=f['root']/'out.tsv'
    write_bam(f['bam'],f['refs'],f['records']+[read(f'edit{i}','C'+'A'*19) for i in range(4)])
    args=r._build_arg_parser().parse_args(['--bam-dir',str(f['root']),
        '--fasta',str(f['fa']),'--config',str(f['flat']),'--bridge-table',str(f['bridge']),
        '--oligo-window-table',str(f['windows']),'--out',str(out)])
    assert r.run(args)==0
    rows=read_rows(out)
    assert len(rows)==1 and rows[0]['label']=='SINGLE_EDIT'
    assert json.loads(rows[0]['allele_counts'])[0]['alt_count']==4


def test_forty_random_haplotypes_have_independent_three_of_five_support(fixture_factory):
    import random
    rng=random.Random(20260920)
    refs,records,cases={},{},[]
    records=[]
    for k in range(40):
        chrom=f'oracle_{k}'
        reference=''.join(rng.choice('ACGT') for _ in range(200))
        refs[chrom]=reference
        kind=k%4
        ref=reference[99:99+(3 if kind==3 else 2 if kind==1 else 1)]
        alt=(''.join(rng.choice([b for b in 'ACGT' if b!=base]) for base in ref) if kind<2
             else ref+'CG' if kind==2 else ref[0])
        allele=r.IntendedAllele(100,ref,alt)
        wt=reference[79:139]
        edited=reference[79:99]+alt+reference[99+len(ref):139]
        cigar='60M' if kind<2 else '21M2I39M' if kind==2 else '21M2D37M'
        for i in range(5):
            records.append(dict(name=f'{chrom}_{i}',contig=chrom,start=79,seq=edited if i<3 else wt,
                                cigar=cigar if i<3 else '60M',flag=0,qual=35,mapq=60))
        cases.append(r.EditSite(chrom,80,139,chrom,(allele,)))
    f=fixture_factory(records=records,extra_refs=refs)
    with pysam.AlignmentFile(str(f['bam'])) as bam,pysam.FastaFile(str(f['fa'])) as fa:
        for es in cases:
            call=r.score_edit_site(bam,fa,es)
            e=call.allele_evidence[0]
            assert (e.alt_count,e.ref_count,e.other_count,e.depth,call.status)==(3,2,0,5,'partial'),es


@pytest.mark.parametrize('allele,seq,cigar', [((305,'AA','A'),'A'*5,'5M1D'),
    ((305,'A','AC'),'A'*5+'C','5M1I')])
def test_terminal_indel_without_observed_right_anchor_is_not_confirmed(fixture_factory,allele,seq,cigar):
    f=fixture_factory(records=[read(f'r{i}',seq,cigar=cigar) for i in range(4)])
    call=score(f,site(allele))
    assert call.status=='nocov' and call.allele_evidence[0].alt_count==0


def test_partial_site_with_an_uncovered_second_allele_remains_low_confidence(fixture_factory):
    f=fixture_factory(records=[read(f'r{i}','C' if i<2 else 'A') for i in range(4)])
    call=score(f,site((301,'A','C'),(310,'A','G')))
    cell=r.reconcile_cell('s',[call],bc_confidence='high')
    assert cell.label=='CONTAMINATION' and cell.confidence=='low'
    assert 'UNCOVERED_ALLELES' in cell.flags


@pytest.mark.parametrize('cigar,seq',[('5M1D1N14M','A'*19),('5M1D1S','A'*6)])
def test_normalized_indel_requires_genomic_not_just_query_right_anchor(fixture_factory,cigar,seq):
    # Independent Stage2 reviewer: adjacent query sequence can be soft-clipped or
    # separated by a reference skip, and cannot corroborate the deletion junction.
    f=fixture_factory(records=[read(f'r{i}',seq,cigar=cigar) for i in range(4)])
    call=score(f,site((305,'AA','A')))
    assert (call.status,call.allele_evidence[0].alt_count,call.cov)==('nocov',0,0)


@pytest.mark.parametrize('mapq,minimum,status,depth,low,unknown',[
    (19,20,'nocov',0,4,0),(20,20,'full',4,0,0),(255,20,'nocov',0,0,4),
    (0,0,'full',4,0,0),(255,0,'nocov',0,0,4)])
def test_variant_mapping_quality_policy_boundaries(fixture_factory,mapq,minimum,status,depth,low,unknown):
    records=[read(f'r{i}','C'+'A'*19) for i in range(4)]
    for record in records:record['mapq']=mapq
    f=fixture_factory(records=records)
    with pysam.AlignmentFile(str(f['bam'])) as bam,pysam.FastaFile(str(f['fa'])) as fa:
        call=r.score_edit_site(bam,fa,site((301,'A','C')),min_mapq=minimum)
    e=call.allele_evidence[0]
    assert (call.status,e.depth,e.alt_count,e.low_mapq_reads,e.unknown_mapq_reads)==(status,depth,depth,low,unknown)


def test_ambiguous_alt_does_not_contaminate_unique_wt_denominator(fixture_factory):
    ambiguous=[read(f'amb{i}','C'+'A'*19) for i in range(4)]
    for record in ambiguous:record['mapq']=0
    f=fixture_factory(records=ambiguous+[read(f'wt{i}','A'*20) for i in range(4)])
    call=score(f,site((301,'A','C')))
    e=call.allele_evidence[0]
    assert (call.status,e.depth,e.ref_count,e.alt_count,e.low_mapq_reads)==('wt',4,4,0,4)
    assert 'AMBIGUOUS_ALIGNMENTS_EXCLUDED' in call.flags


@pytest.mark.parametrize('allele,seq,cigar', [
    ((305,'A','AC'),'AAAAGC'+'A'*15,'5M1I15M'),
    ((305,'AA','A'),'AAAAG'+'A'*14,'5M1D14M'),
])
def test_normalized_indel_does_not_replace_observed_mutant_anchor_with_reference(fixture_factory,allele,seq,cigar):
    # Declared ALT starts with A. The aligned allele starts with G; FASTA cannot
    # substitute A for that observed G when constructing the candidate CIGAR event.
    f=fixture_factory(records=[read(f'r{i}',seq,cigar=cigar) for i in range(4)])
    call=score(f,site(allele))
    e=call.allele_evidence[0]
    assert (call.status,e.alt_count,e.ref_count,e.other_count,e.depth)==('other',0,0,4,4)


@pytest.mark.parametrize('allele,seq,cigar,expected',[
    ((305,'AA','A'),'A'*5+'G'+'A'*13,'8M1D11M','other'),
    ((305,'AA','A'),'A'*6+'G'+'A'*12,'8M1D11M','other'),
    ((305,'A','AA'),'A'*5+'G'+'A'*15,'8M1I12M','other'),
])
def test_shifted_indel_cannot_cross_an_observed_snp_interrupting_the_repeat(fixture_factory,allele,seq,cigar,expected):
    # In FASTA this entire interval is A. In the observed molecule G interrupts
    # the run between the requested event and the CIGAR indel. Moving the indel
    # across G changes the haplotype, so reference-only normalization is invalid.
    f=fixture_factory(records=[read(f'r{i}',seq,cigar=cigar) for i in range(4)])
    call=score(f,site(allele))
    e=call.allele_evidence[0]
    assert call.status==expected and e.alt_count==0 and e.depth==4
    assert (e.other_count,e.ref_count)==((4,0) if expected=='other' else (0,4))


@pytest.mark.parametrize('allele,seq,cigar',[
    ((305,'A','AA'),'A'*22,'5M1I3M1I12M'),
    ((305,'AA','A'),'A'*18,'5M1D3M1D10M'),
])
def test_two_equivalent_indels_are_not_one_exact_intended_haplotype(fixture_factory,allele,seq,cigar):
    # Each individual event left-normalizes to the target, but the molecule has
    # two changes. Its spanning haplotype cannot equal a design containing one.
    f=fixture_factory(records=[read(f'r{i}',seq,cigar=cigar) for i in range(4)])
    e=score(f,site(allele)).allele_evidence[0]
    assert (e.alt_count,e.other_count,e.depth)==(0,4,4)


@pytest.mark.parametrize('cigar',[ '8M1D1P11M','8M1P1D11M'])
def test_cigar_padding_does_not_remove_real_genomic_indel_anchors(fixture_factory,cigar):
    # P consumes neither query nor reference sequence; aligned bases still flank
    # this repeat-equivalent event on both genomic sides.
    f=fixture_factory(records=[read(f'r{i}','A'*19,cigar=cigar) for i in range(4)])
    call=score(f,site((305,'AA','A')))
    assert (call.status,call.allele_evidence[0].alt_count,call.cov)==('full',4,4)


@pytest.mark.parametrize('allele',[(302,'A','C'),(301,'AA','AC'),(302,'AA','CA'),(301,'AAA','ACA')])
def test_unchanged_ref_alt_padding_does_not_change_variant_coverage(fixture_factory,allele):
    # All four declarations specify the same A>C at302. Only that changing base
    # is covered; unchanged representational padding is not missing edit evidence.
    f=fixture_factory(records=[read(f'r{i}','C',start=301) for i in range(4)])
    call=score(f,site(allele))
    e=call.allele_evidence[0]
    assert (call.status,e.alt_count,e.ref_count,e.other_count,e.depth)==('full',4,0,0,4)
    assert (e.allele.pos,e.allele.ref,e.allele.alt)==allele




def test_two_independent_designed_sites_are_multi_edit(fixture_factory):
    records=[read(f'left{i}','C',start=300) for i in range(4)]+[
        read(f'right{i}','G',start=400) for i in range(4)]
    f=fixture_factory(records=records)
    left=score(f,site((301,'A','C')))
    right=score(f,r.EditSite('chrII',401,410,'oligo_2',(r.IntendedAllele(401,'A','G'),)))
    cell=r.reconcile_cell('s',[left,right],bc_confidence='high')
    assert (cell.label,cell.n_full,cell.n_partial)==('MULTI_EDIT',2,0)
    assert [(s.allele_evidence[0].alt_count,s.cov) for s in cell.sites]==[(4,4),(4,4)]


@pytest.mark.parametrize('allele,seq,start', [((301,'AA','CG'),'C',300),((301,'A','AC'),'A',300)])
def test_incomplete_required_variant_context_is_unknown(fixture_factory,allele,seq,start):
    f=fixture_factory(records=[read(f'r{i}',seq,start=start) for i in range(4)])
    call=score(f,site(allele))
    assert (call.status,call.cov)==('nocov',0)


@pytest.mark.parametrize('bad', ['ambiguous','low_quality'])
def test_insertion_base_itself_must_be_called_and_high_quality(fixture_factory,bad):
    f=fixture_factory(records=[read('r','A'*5+('N' if bad=='ambiguous' else 'C')+'A'*15,cigar='5M1I15M')])
    with pysam.AlignmentFile(str(f['bam'])) as bam:
        molecule=next(bam.fetch())
    if bad=='low_quality':
        qualities=list(molecule.query_qualities);qualities[5]=0;molecule.query_qualities=qualities
    assert r._read_allele(molecule,r.IntendedAllele(305,'A','AC'),13) is None


def test_absent_sequence_or_qualities_cannot_supply_an_allele():
    for seq,qualities in [(None,None),('A'*20,None)]:
        molecule=SimpleNamespace(reference_start=300,reference_end=320,query_sequence=seq,query_qualities=qualities,
                                 cigartuples=[(0,20)])
        assert r._read_allele(molecule,r.IntendedAllele(301,'A','C'),13) is None
        assert r._equivalent_indel(molecule,None,'chrII',r.IntendedAllele(305,'AA','A'),(305,'AA','A'),13,{}) is None


def test_unsupported_cigar_does_not_invent_a_sequence():
    molecule=SimpleNamespace(reference_start=300,reference_end=320,query_sequence='A'*20,query_qualities=[35]*20,
                             cigartuples=[(0,5),(9,1),(0,15)])
    assert r._read_allele(molecule,r.IntendedAllele(306,'A','C'),13) is None


@pytest.mark.parametrize('prefix,pos,expected', [('',1,(1,'A','AA')),('N',5,(2,'A','AA'))])
def test_indel_left_normalization_stops_at_contig_start_or_unknown_reference(fixture_factory,prefix,pos,expected):
    f=fixture_factory(records=[],extra_refs={'boundary':prefix+'A'*100})
    with pysam.FastaFile(str(f['fa'])) as fa:
        assert r._normalize_allele(fa,'boundary',r.IntendedAllele(pos,'A','AA'))==expected


def test_shifted_equivalence_with_unobserved_intervening_context_is_unknown(fixture_factory):
    # N crosses the context between requested deletion and observed repeat-shifted
    # deletion. A downstream real junction cannot restore the missing phase.
    f=fixture_factory(records=[read(f'r{i}','A'*19,cigar='5M1N3M1D11M') for i in range(4)])
    call=score(f,site((305,'AA','A')))
    assert (call.status,call.allele_evidence[0].alt_count,call.cov)==('nocov',0,0)


@pytest.mark.parametrize('kwargs',[{'min_mapq':-1},{'min_mapq':255},{'min_mapq':True},{'min_mapq':1.5},
    {'min_base_qual':-1},{'min_base_qual':1.5}])
def test_invalid_variant_filter_policies_are_rejected(fixture_factory,kwargs):
    f=fixture_factory(records=[])
    with pysam.AlignmentFile(str(f['bam'])) as bam,pysam.FastaFile(str(f['fa'])) as fa:
        with pytest.raises(ValueError):r.score_edit_site(bam,fa,site((301,'A','C')),**kwargs)


@pytest.mark.parametrize('text', ['', 'oligo_name\tchrom\tdonor_start_coord\tdonor_end_coord\n',
    'oligo_name\tchrom\tchrom\tdonor_start_coord\tdonor_end_coord\na\tc\tc\t1\t2\n',
    'name\tchrom\tstart\tend\na\tc\t1\t2\n'])
def test_missing_duplicate_and_empty_window_headers_fail(tmp_path,text):
    p=tmp_path/'windows.tsv';p.write_text(text)
    with pytest.raises(ValueError):r.load_oligo_windows(str(p))


def test_unknown_site_status_is_not_silently_interpreted():
    with pytest.raises(ValueError,match='unknown site status'):
        r.reconcile_cell('s',[r.SiteCall(site(),1,4,'unexpected')])


def test_build_sites_rejects_reference_mismatch_before_deduplication(fixture_factory):
    f=fixture_factory(records=[])
    with pysam.FastaFile(str(f['fa'])) as fa:
        with pytest.raises(ValueError,match='REF does not match'):
            r.build_edit_sites(['bad'],{'bad':site((301,'G','C'))},fasta=fa)


def test_different_contigs_are_independent_and_conflicting_window_scopes_fail():
    a=site((301,'A','C'))
    b=r.EditSite('other',301,320,'b',(r.IntendedAllele(301,'A','C'),))
    assert len(r.build_edit_sites(['a','b'],{'a':a,'b':b}))==2
    with pytest.raises(ValueError,match='conflicting intended designs'):
        r.build_edit_sites(['a','legacy'],{'a':a,'legacy':site()})


def test_cli_entrypoint_and_empty_suffix(fixture_factory):
    f=fixture_factory();out=f['root']/'cli.tsv'
    args=['--bam-dir',str(f['root']),'--fasta',str(f['fa']),'--config',str(f['flat']),
          '--bridge-table',str(f['bridge']),'--oligo-window-table',str(f['windows']),'--out',str(out)]
    with pytest.raises(ValueError,match='nonempty'):r.main(args+['--bam-suffix',''])
    assert not out.exists()
    assert r.main(args)==0
    rows=read_rows(out)
    assert len(rows)==1 and rows[0]['label']=='LOW_CONF'
@pytest.mark.parametrize('allele,seq,cigar', [
    ((305,'AAA','AC'),'A'*5+'C'+'A'*14,'5M2D1I14M'),
    ((305,'AAA','AC'),'A'*5+'C'+'A'*14,'5M1I2D14M'),
    ((305,'AA','ACG'),'A'*5+'CG'+'A'*14,'5M1D2I14M'),
    ((304,'AAAA','AAC'),'A'*5+'C'+'A'*14,'5M2D1I14M'),
])
def test_unchanged_padding_trim_preserves_required_complex_indel_anchor(fixture_factory,allele,seq,cigar):
    # Independently specified anchored delins: deleting AA and inserting C (or
    # deleting A and inserting CG) must survive trimming of unchanged padding.
    f=fixture_factory(records=[read(f'r{i}',seq,cigar=cigar) for i in range(4)])
    call=score(f,site(allele))
    e=call.allele_evidence[0]
    assert (call.status,e.alt_count,e.ref_count,e.other_count,e.depth)==('full',4,0,0,4)
    assert (e.allele.pos,e.allele.ref,e.allele.alt)==allele
@pytest.mark.parametrize('kind,start,end,expected_vaf,expected_cov', [
    ('low_mapq',301,301,0,4),
    ('unknown_mapq',301,301,0,4),
    ('refskip',306,306,0,0),
    ('deletion',306,306,1,4),
    ('ambiguous',301,301,0,0),
])
def test_legacy_window_observables_keep_missing_and_rejected_evidence_explicit(fixture_factory,kind,start,end,expected_vaf,expected_cov):
    if kind in ('low_mapq','unknown_mapq'):
        records=[read(f'wt{i}','A') for i in range(4)]
        records += [dict(read(f'bad{i}','C'),mapq=19 if kind=='low_mapq' else 255) for i in range(4)]
    elif kind in ('refskip','deletion'):
        records=[read(f'r{i}','A'*19,cigar='5M1'+('N' if kind=='refskip' else 'D')+'14M') for i in range(4)]
    else:
        records=[read(f'r{i}','N') for i in range(4)]
    f=fixture_factory(records=records)
    call=score(f,r.EditSite('chrII',start,end,'legacy'))
    assert (call.vaf,call.cov,call.status)==(expected_vaf,expected_cov,'unassessed')
    assert call.allele_evidence==()
    assert r.reconcile_cell('s',[call],bc_confidence='high').label=='LOW_CONF'


def test_conflicting_rows_for_one_oligo_are_rejected(tmp_path):
    p=tmp_path/'windows.tsv'
    p.write_text('oligo_name\tchrom\tdonor_start_coord\tdonor_end_coord\n'
                 'a\tchrII\t301\t310\na\tchrII\t301\t311\n')
    with pytest.raises(ValueError,match='conflicting donor windows'):
        r.load_oligo_windows(str(p))
