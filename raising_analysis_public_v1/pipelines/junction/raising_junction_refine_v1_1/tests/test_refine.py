"""Synthetic tests only; no claim of real sample accuracy."""
import csv, gzip, hashlib, json, random, subprocess, sys, tempfile, unittest
from pathlib import Path
PACKAGE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(PACKAGE))
sys.path.insert(0,str(PACKAGE/'support'))
import junction_core as exact
import refine_templates as templates
import refine_alignment as alignment
import raising_junction_refine as runner
import pilot_io

def rc(s):
    return s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]

def random_dna(n, seed=183):
    rng=random.Random(seed)
    return ''.join(rng.choice('ACGT') for _ in range(n))

def event(identifier='J0',boundary=500,strand='+',deletion=10,insertion='',support=10,side='3p',mh=0,chrom='chrTest'):
    return dict(sample='Synthetic_3prm_S1',library_side=side,event_id=identifier,chrom=chrom,host_strand=strand,
      host_boundary0=boundary,host_boundary_min0=boundary-mh if strand=='+' else boundary,
      host_boundary_max0=boundary if strand=='+' else boundary+mh,
      donor_cut_oriented_min0=800-deletion-mh,donor_cut_oriented_max0=800-deletion,
      donor_deletion_min_bp=deletion,donor_deletion_max_bp=deletion+mh,mh_length_bp=mh,
      mh_sequence=(DONOR if side=='3p' else rc(DONOR))[800-deletion-mh:800-deletion],
      insertion_length_bp=len(insertion),insertion_sequence=insertion,pair_support=support,
      simple_pair_support=support,prior_dd_pair_support=0,complex_host_pair_support=0)

def make_templates(rows=None,side='3p',genome=None):
    rows=[event(side=side)] if rows is None else rows
    genome=random_dna(2000) if genome is None else genome
    refs={r['chrom']:genome for r in rows}
    ts=templates.build_templates(templates.cluster_candidates(rows),DONOR,side,
        lambda c,s,e:refs[c][s:e],lambda c:len(refs[c]))
    return ts,refs

def call(t,query=None,quality=None,intervals=None):
    query=t['sequence'] if query is None else query
    quality='I'*len(query) if quality is None else quality
    intervals=[(t['host_interval_start0'],t['host_interval_end0'])] if intervals is None else intervals
    return alignment.align_to_template(query,quality,t,DONOR,intervals)

class TemplateTests(unittest.TestCase):
    def test_donor_length_and_primer_positions(self):
        self.assertEqual(len(DONOR),800)
        self.assertEqual(exact.primer_info(DONOR,'3p')['outward_bases'],59)
        self.assertEqual(exact.primer_info(DONOR,'5p')['outward_bases'],65)

    def test_cluster_span_does_not_chain(self):
        cs=templates.cluster_candidates([event('a',100),event('b',108),event('c',116),event('d',100,strand='-'),event('e',100,chrom='MT')],max_span=10)
        plus=[c for c in cs if c['chrom']=='chrTest' and c['host_strand']=='+']
        self.assertEqual([sorted(r['event_id'] for r in c['rows']) for c in plus],[['a','b'],['c']])
        self.assertTrue(all(c['boundary_max0']-c['boundary_min0']<=10 for c in cs))

    def test_plus_has_common_human_endpoint(self):
        ts,refs=make_templates([event('a',500),event('b',505)])
        self.assertEqual({t['host_interval_end0'] for t in ts},{625})
        for t in ts:
            self.assertEqual(t['sequence'][t['host_start_t0']:],refs['chrTest'][t['host_interval_start0']:625])

    def test_minus_has_common_reverse_complement_interval(self):
        ts,refs=make_templates([event('a',500,strand='-'),event('b',505,strand='-')])
        self.assertEqual({t['host_interval_start0'] for t in ts},{380})
        for t in ts:
            self.assertEqual(t['sequence'][t['host_start_t0']:],rc(refs['chrTest'][380:t['host_interval_end0']]))

    def test_five_prime_outward_orientation(self):
        ts,_=make_templates([event(side='5p',deletion=12)],side='5p')
        t=ts[0]
        self.assertEqual(t['sequence'][:t['donor_end_t0']],rc(DONOR)[735:788])

    def test_sequence_equivalence_keeps_true_delta_one_distinct(self):
        ts,_=make_templates([event('intact',deletion=0),event('equivalent',deletion=1,insertion=DONOR[-1]),event('delta1',deletion=1)])
        by={t['template_id']:t for t in ts}
        self.assertEqual(len(ts),3)
        self.assertEqual(by['intact']['equivalent_template_group'],by['equivalent']['equivalent_template_group'])
        self.assertNotEqual(by['intact']['equivalent_template_group'],by['delta1']['equivalent_template_group'])

    def test_real_donor_tail_snv_is_alternative_not_forced_merge(self):
        rows=[event('major',deletion=0,support=220),event('minor',deletion=10,insertion='ACTGTCGAAT',support=1)]
        review=templates.build_donor_only_review(rows,DONOR,'3p')
        minor=next(r for r in review if r.get('template_id',r.get('event_id'))=='minor')
        self.assertEqual(minor['relationship'],'substitution_only_alternative')
        self.assertTrue(minor['limited_host_context'])

class AlignmentTests(unittest.TestCase):
    def test_exact_supported(self):
        r=call(make_templates()[0][0])
        self.assertEqual(r['status'],'supported',r)

    def test_internal_snv_allowed(self):
        t=make_templates()[0][0];s=t['sequence']
        q=s[:15]+next(c for c in 'ACGT' if c!=s[15])+s[16:]
        r=call(t,q)
        self.assertEqual(r['status'],'supported',r)

    def test_internal_small_insertion_allowed(self):
        t=make_templates()[0][0];s=t['sequence']
        r=call(t,s[:15]+'AA'+s[15:])
        self.assertEqual(r['status'],'supported',r)

    def test_boundary_snv_review(self):
        t=make_templates()[0][0];s=t['sequence'];p=t['donor_end_t0']-1
        q=s[:p]+next(c for c in 'ACGT' if c!=s[p])+s[p+1:]
        r=call(t,q)
        self.assertEqual(r['status'],'review',r)

    def test_low_quality_junction_rejected(self):
        t=make_templates()[0][0];q=list('I'*len(t['sequence']));q[t['donor_end_t0']]='!'
        r=call(t,quality=''.join(q))
        self.assertEqual(r['status'],'rejected',r)

    def test_mate_only_not_junction(self):
        t=make_templates()[0][0]
        for q in (t['sequence'][:t['donor_end_t0']],t['sequence'][t['host_start_t0']:]):
            self.assertEqual(call(t,q)['status'],'rejected')

    def test_wrong_global_locus_rejected(self):
        self.assertEqual(call(make_templates()[0][0],intervals=[(1200,1300)])['status'],'rejected')

    def test_reverse_host_and_five_prime_alignment_supported(self):
        for side in ('3p','5p'):
            t=make_templates([event(side=side,strand='-')],side=side)[0][0]
            self.assertEqual(call(t)['status'],'supported')

    def test_mh_touching_variant_keeps_reference_mh_separate(self):
        genome=random_dna(2000)
        genome=genome[:495]+DONOR[785:790]+genome[500:]
        t=make_templates([event(mh=5)],genome=genome)[0][0]
        s=t['sequence'];p=t['donor_end_t0']-4
        q=s[:p]+next(c for c in 'ACGT' if c!=s[p])+s[p+1:]
        r=call(t,q)
        self.assertEqual(r['status'],'review',r)
        self.assertEqual(r['best_alignment']['template_reference_mh_length_bp'],5)
        self.assertNotIn('observed_mh_length_bp',r['best_alignment'])

    def test_minor_snv_allele_fits_itself_better_than_abundant_template(self):
        ts,_=make_templates([event('major',deletion=0,support=220),event('minor',deletion=10,insertion='ACTGTCGAAT',support=1)])
        by={t['template_id']:t for t in ts};q=by['minor']['sequence']
        own=call(by['minor'],q);major=call(by['major'],q)
        self.assertEqual(own['status'],'supported',own)
        self.assertEqual(major['status'],'supported',major)
        self.assertGreater(own['score'],major['score'])
        self.assertEqual(major['best_alignment']['substitution_bases'],1)

    def test_repeated_boundary_gap_placements_remain_ambiguous(self):
        genome=random_dna(2000)
        genome=genome[:500]+'AAAAAA'+genome[506:]
        t=make_templates([event(deletion=0)],genome=genome)[0][0]
        s=t['sequence'];p=t['host_start_t0']
        r=call(t,s[:p]+s[p+1:])
        self.assertEqual(r['status'],'ambiguous_alignment',r)

    def test_true_terminal_deletion_prefers_own_template(self):
        ts,_=make_templates([event('major',deletion=0,support=220),event('delta1',deletion=1,support=1)])
        by={t['template_id']:t for t in ts};q=by['delta1']['sequence']
        own=call(by['delta1'],q);other=call(by['major'],q)
        self.assertEqual(own['status'],'supported',own)
        self.assertIn(other['status'],('review','rejected','ambiguous_alignment'))
        self.assertGreater(own['score'],other['score'])


def write_table(path,rows):
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'wt',newline='') as f:
        w=csv.DictWriter(f,list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)

def read_table(path):
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt',newline='') as f:
        return list(csv.DictReader(f,delimiter='\t'))

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def write_fasta(path,records):
    entries=[]
    with path.open('wb') as f:
        for name,s in records.items():
            f.write(('>'+name+'\n').encode());offset=f.tell()
            f.write((s+'\n').encode())
            entries.append(f'{name}\t{len(s)}\t{offset}\t{len(s)}\t{len(s)+1}\n')
    Path(str(path)+'.fai').write_text(''.join(entries))

def make_runner_fixture(root):
    sample='Synthetic_3prm_S1'
    evidence=root/'results'/(sample+'_fourway_full_v1')
    pilot=root/'results'/(sample+'_junction_pilot_v1')
    raw=root/'vendor_results'/'Raw_fastq'
    ref=root/'ref'/'synthetic.fa'
    for p in (evidence,pilot,raw,ref.parent):p.mkdir(parents=True,exist_ok=True)
    genomes={'chrTest':random_dna(2000),'chrOther':random_dna(2000,331),'MT':random_dna(2000,777), '__RAISING_DONOR_800__':DONOR}
    write_fasta(ref,genomes)
    dp=DONOR[741:790];hp=genomes['chrTest'][500:570]
    mhp=hp[:25]+next(c for c in 'ACGT' if c!=hp[25])+hp[26:]
    mt=genomes['MT'][500:570]
    queries={
      'pairA':{'R1':dp+hp,'R2':rc(dp+hp)},
      'pairB':{'R1':dp+mhp,'R2':hp},
      'pairC':{'R1':dp,'R2':hp},
      'pairD':{'R1':dp,'R2':genomes['chrOther'][500:570]},
      'pairE':{'R1':dp,'R2':rc(dp)},
      'pairM':{'R1':dp+mt,'R2':mt}}
    raw_records={};paths={}
    for mate in ('R1','R2'):
        path=raw/f'{sample}_L001_{mate}_001.fastq.gz';paths[mate]=path
        with gzip.open(path,'wt') as f:
            for rid,ends in queries.items():
                s='CGT'+ends[mate]+'AC';record=(f'@{rid}/{mate[-1]}',s,'+','I'*len(s))
                raw_records[(rid,mate)]=record;f.write('\n'.join(record)+'\n')
    pairs=[];reads=[];hits=[]
    for rid,ends in queries.items():
        pairs.append(dict(read_id=rid,classification='other_unresolved' if rid=='pairE' else 'host_mappable_other',host_mapping_supported=int(rid!='pairE'),donor_rearrangement_supported=0))
        for mate,q in ends.items():
            reads.append(dict(read_id=rid,mate=mate,raw_length=len(q)+5,cleaned_length=len(q),cleaned_raw_start_1based=4,cleaned_raw_end_1based=len(q)+3))
    for rid,mate,chrom,cigar,strand in (
      ('pairA','R1','chrTest','49S70M','+'),('pairA','R2','chrTest','70M49S','-'),
      ('pairB','R1','chrTest','49S70M','+'),('pairB','R2','chrTest','70M','+'),
      ('pairC','R2','chrTest','70M','+'),('pairD','R2','chrOther','70M','+'),
      ('pairM','R1','MT','49S70M','+'),('pairM','R2','MT','70M','+')):
        hits.append(dict(read_id=rid,mate=mate,reference=chrom,position_1based=501,cigar=cigar,strand=strand,accepted=1,mapq=60,supplementary=0))
    for name,rows in [('pair_classification',pairs),('read_evidence',reads),('host_alignment_evidence',hits)]:
        write_table(evidence/(name+'.tsv.gz'),rows)
    digest=hashlib.sha256()
    for rid in queries:pilot_io.update_stream_hash(digest,rid,raw_records[(rid,'R1')],raw_records[(rid,'R2')])
    manifest=dict(sample=sample,side='3p',status='complete',total_pairs=6,host_supported_pairs=5,analysis_version='1.0.0',
       inputs=dict(lanes=[{mate:dict(path=str(p),sha256=sha(p)) for mate,p in paths.items()}],normalized_paired_stream_sha256=digest.hexdigest()),
       reference=dict(combined_fasta=str(ref),signature=dict(donor_sha256=hashlib.sha256(DONOR.encode()).hexdigest(),donor_contig='__RAISING_DONOR_800__')),
       final_output_sha256={p.name:sha(p) for p in evidence.glob('*.tsv.gz')})
    manifest_path=evidence/'run_summary.json';manifest_path.write_text(json.dumps(manifest))
    jr=[event('J0',support=1),event('Jmt',support=1,chrom='MT')]
    for row in jr:row['example_candidate_json']='{}'
    write_table(pilot/'junctions.tsv',jr)
    old=[]
    for rid in ['pairA','pairB','pairC','pairD','pairM']:
        good=rid in ('pairA','pairM')
        diagnostic=[dict(diagnostics=dict(rejected_counts={'no_exact_host_anchor':2})),dict(diagnostics=dict(rejected_counts={'no_exact_host_anchor':1}))]
        old.append(dict(read_id=rid,status='called_simple_candidate' if good else 'unresolved',event_id='J0' if rid=='pairA' else 'Jmt' if rid=='pairM' else '',diagnostics_json=json.dumps([] if good else diagnostic),complex_host_evidence=0,prior_donor_rearrangement_supported=0))
    write_table(pilot/'pair_junction_status.tsv.gz',old)
    ps=dict(analysis_type='raising_exact_anchor_junction_pilot',analysis_version='0.1.0',sample=sample,library_side='3p',status='complete',validation={'all':True},
       stats=dict(observed_junction_patterns_total=2,host_supported_pairs=5,called_pairs_total=2),
       output_sha256={p.name:sha(p) for p in pilot.iterdir() if p.is_file()},
       upstream_input_sha256={p.name:sha(p) for p in evidence.iterdir() if p.is_file()},
       reference=dict(path=str(ref),sha256=sha(ref),donor_sha256=hashlib.sha256(DONOR.encode()).hexdigest()))
    (pilot/'summary.json').write_text(json.dumps(ps))
    return sample,evidence,pilot,paths


class RunnerTests(unittest.TestCase):
    def invoke(self,root,sample,*extra):
        return subprocess.run([sys.executable,str(PACKAGE/'raising_junction_refine.py'),'--root',str(root),'--sample',sample,*extra],capture_output=True,text=True,timeout=30)

    def test_reasons_count_pairs_not_orientation_attempts(self):
        attempt={'diagnostics':{'rejected_counts':{'reasonA':3,'reasonB':1}}}
        rows=[{'status':'unresolved','diagnostics_json':json.dumps([attempt,attempt])},
              {'status':'unresolved','diagnostics_json':json.dumps([{'diagnostics':{'rejected_counts':{'reasonA':1}}}])},
              {'status':'ambiguous','diagnostics_json':'[]'}]
        result=runner.summarize_upstream_reasons(rows)
        by={(r['old_status'],r['reason']):r for r in result}
        self.assertEqual(by[('unresolved','reasonA')]['affected_pairs'],2)
        self.assertEqual(by[('unresolved','reasonB')]['affected_pairs'],1)
        self.assertEqual(by[('unresolved','reasonB')]['denominator_pairs'],2)
        self.assertTrue(all(r['nonexclusive']==1 for r in result))

    def result_record(self,mate,group,score=100,status='supported',cluster='c1',context='ctx1'):
        return dict(mate=mate,cluster_id=cluster,comparison_context_id=context,equivalent_template_group=group,template_id=group,
                    result=dict(status=status,score=score,best_alignment={'flags':[]}))

    def test_mate_intersection_can_resolve_template_tie(self):
        rows=[self.result_record('R1','a'),self.result_record('R1','b'),self.result_record('R2','a')]
        d=runner.resolve_pair(rows,{'a':['a'],'b':['b']})
        self.assertEqual(d['status'],'template_supported')
        self.assertEqual(d['template_ids'],['a'])
        self.assertEqual(d['supporting_mates'],['R1','R2'])

    def test_discordant_mates_remain_ambiguous(self):
        rows=[self.result_record('R1','a'),self.result_record('R2','b')]
        self.assertEqual(runner.resolve_pair(rows,{'a':['a'],'b':['b']})['status'],'template_ambiguous')

    def test_cross_locus_qualifying_hits_not_reduced_to_global_best(self):
        rows=[self.result_record('R1','a',score=200,cluster='c1',context='ctx1'),self.result_record('R1','b',score=100,cluster='c2',context='ctx2')]
        self.assertEqual(runner.resolve_pair(rows,{'a':['a'],'b':['b']})['status'],'template_ambiguous')

    def test_unique_template_with_boundary_variant_stays_review(self):
        rows=[self.result_record('R1','a',status='review'),self.result_record('R2','a')]
        self.assertEqual(runner.resolve_pair(rows,{'a':['a']})['status'],'junction_variant_review')

    def test_display_clusters_are_not_separate_competitions(self):
        rows=[self.result_record('R1','major',score=200,cluster='display_left',context='same_hit'),
              self.result_record('R1','minor',score=194,cluster='display_right',context='same_hit')]
        d=runner.resolve_pair(rows,{'major':['major'],'minor':['minor']})
        self.assertEqual(d['status'],'template_supported')
        self.assertEqual(d['template_ids'],['major'])

    def test_distinct_accepted_contexts_preserved_even_if_display_id_shared(self):
        rows=[self.result_record('R1','major',score=200,cluster='display',context='chr1_hit'),
              self.result_record('R1','minor',score=194,cluster='display',context='chr2_hit')]
        self.assertEqual(runner.resolve_pair(rows,{'major':['major'],'minor':['minor']})['status'],'template_ambiguous')

    def test_accepted_contexts_merge_only_overlapping_same_reference_orientation(self):
        hits=[{'reference':'chr1','strand':'+','intervals0':[(100,160)]},
              {'reference':'chr1','strand':'+','intervals0':[(150,190)]},
              {'reference':'chr1','strand':'+','intervals0':[(300,350)]},
              {'reference':'chr2','strand':'+','intervals0':[(100,160)]},
              {'reference':'chr1','strand':'-','intervals0':[(100,160)]}]
        contexts=runner.accepted_host_contexts(hits)
        self.assertEqual(len(contexts),4)
        plus=sorted((c['span_start0'],c['span_end0']) for c in contexts if c['chrom']=='chr1' and c['host_strand']=='+')
        self.assertEqual(plus,[(100,190),(300,350)])
        flipped=runner.accepted_host_contexts(hits,flip=True)
        minus=sorted((c['span_start0'],c['span_end0']) for c in flipped if c['chrom']=='chr1' and c['host_strand']=='-')
        self.assertEqual(minus,plus)

    def test_shifted_snv_minor_does_not_gain_longer_template_advantage(self):
        genome=random_dna(2000,183)
        fetch=lambda chrom,a,b:genome[a:b]
        for strand in ('+','-'):
            for shift in (12,18,19,20):
                with self.subTest(strand=strand,shift=shift):
                    human=genome[500:700] if strand=='+' else rc(genome[300:500])
                    intervening=human[:shift]
                    intervening=intervening[:5]+next(c for c in 'ACGT' if c!=intervening[5])+intervening[6:]
                    rows=[event('major',boundary=500,strand=strand,deletion=0,support=220),
                          event('minor',boundary=500+shift if strand=='+' else 500-shift,strand=strand,deletion=0,insertion=intervening,support=1)]
                    query=DONOR[741:]+human
                    decisions=[]
                    for span in (0,10,20):
                        display=templates.build_templates(templates.cluster_candidates(rows,max_span=span),DONOR,'3p',fetch,lambda chrom:len(genome))
                        runner.assign_sequence_equivalence(display,DONOR,'3p',fetch,lambda chrom:len(genome))
                        by_id={t['template_id']:t for t in display}
                        intervals=[(500,700)] if strand=='+' else [(300,500)]
                        context=runner.accepted_host_contexts([{'reference':'chrTest','strand':strand,'intervals0':intervals}])[0]
                        common=runner.build_comparison_templates(context,runner.template_index(display),by_id,'major',DONOR,'3p',fetch,lambda chrom:len(genome))
                        self.assertEqual({t['template_id'] for t in common},{'major','minor'})
                        endpoints={t['host_interval_end0'] if strand=='+' else t['host_interval_start0'] for t in common}
                        self.assertEqual(len(endpoints),1)
                        scores={};records=[];members={}
                        for t in common:
                            result=alignment.align_to_template(query,'I'*len(query),t,DONOR,intervals)
                            scores[t['template_id']]=result['score']
                            self.assertNotEqual(result['status'],'rejected',result)
                            members.setdefault(t['equivalent_template_group'],[]).append(t['template_id'])
                            records.append({'mate':'R1','cluster_id':t['cluster_id'],
                                'comparison_context_id':t['comparison_context_id'],
                                'equivalent_template_group':t['equivalent_template_group'],
                                'template_id':t['template_id'],'result':result})
                        self.assertEqual(scores['major']-scores['minor'],6)
                        decision=runner.resolve_pair(records,members)
                        self.assertEqual(decision['status'],'template_supported',decision)
                        self.assertEqual(decision['template_ids'],['major'],decision)
                        decisions.append((decision['status'],decision['template_ids']))
                    self.assertTrue(all(d==decisions[0] for d in decisions))

    def test_equal_representations_do_not_change_decision_with_display_span(self):
        genome=random_dna(2000)
        fetch=lambda chrom,a,b:genome[a:b]
        rows=[event('major',deletion=0,support=220),
              event('equivalent_shift',boundary=512,deletion=0,insertion=genome[500:512],support=1)]
        query=DONOR[741:]+genome[500:700]
        decisions=[]
        for span in (0,10,20):
            display=templates.build_templates(templates.cluster_candidates(rows,max_span=span),DONOR,'3p',fetch,lambda c:len(genome))
            runner.assign_sequence_equivalence(display,DONOR,'3p',fetch,lambda c:len(genome))
            by={t['template_id']:t for t in display}
            context=runner.accepted_host_contexts([{'reference':'chrTest','strand':'+','intervals0':[(500,700)]}])[0]
            common=runner.build_comparison_templates(context,runner.template_index(display),by,'major',DONOR,'3p',fetch,lambda c:len(genome))
            self.assertEqual(len({t['sequence'] for t in common}),1)
            records=[];members={}
            for t in common:
                result=alignment.align_to_template(query,'I'*len(query),t,DONOR,[(500,700)])
                members.setdefault(t['equivalent_template_group'],[]).append(t['template_id'])
                records.append(dict(mate='R1',template_id=t['template_id'],cluster_id=t['cluster_id'],
                    comparison_context_id=t['comparison_context_id'],equivalent_template_group=t['equivalent_template_group'],result=result))
            d=runner.resolve_pair(records,members)
            decisions.append((d['status'],sorted(d['template_ids'])))
        self.assertTrue(all(d==decisions[0] for d in decisions),decisions)

    def test_unequal_comparison_endpoints_are_rejected(self):
        a=self.result_record('R1','a');b=self.result_record('R1','b')
        a['comparison_common_host_endpoint0']=620;b['comparison_common_host_endpoint0']=632
        with self.assertRaises(ValueError):
            runner.resolve_pair([a,b],{'a':['a'],'b':['b']})

    def test_equal_distinct_template_scores_remain_ambiguous(self):
        rows=[self.result_record('R1','a'),self.result_record('R1','b')]
        self.assertEqual(runner.resolve_pair(rows,{'a':['a'],'b':['b']})['status'],'template_ambiguous')

    def test_minus_equivalent_alleles_share_global_identity(self):
        genome=random_dna(2000)
        fetch=lambda chrom,a,b:genome[a:b]
        rows=[event('major',boundary=500,strand='-',deletion=0),
              event('equivalent',boundary=488,strand='-',deletion=0,insertion=rc(genome[488:500]))]
        ts=templates.build_templates(templates.cluster_candidates(rows,max_span=0),DONOR,'3p',fetch,lambda c:len(genome))
        runner.assign_sequence_equivalence(ts,DONOR,'3p',fetch,lambda c:len(genome))
        self.assertEqual(len({t['equivalent_template_group'] for t in ts}),1)

    def test_global_identity_does_not_erase_true_terminal_one_base_deletion(self):
        genome=random_dna(2000)
        fetch=lambda chrom,a,b:genome[a:b]
        for side in ('3p','5p'):
            for strand in ('+','-'):
                rows=[event('major',strand=strand,side=side,deletion=0),event('delta1',strand=strand,side=side,deletion=1)]
                ts=templates.build_templates(templates.cluster_candidates(rows,max_span=0),DONOR,side,fetch,lambda c:len(genome))
                runner.assign_sequence_equivalence(ts,DONOR,side,fetch,lambda c:len(genome))
                self.assertEqual(len({t['equivalent_template_group'] for t in ts}),2,(side,strand))

    def test_end_to_end_rescue_dedup_mt_and_unresolved(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);sample,evidence,pilot,paths=make_runner_fixture(root)
            inputs=[p for folder in (evidence,pilot,root/'ref',root/'vendor_results') for p in folder.rglob('*') if p.is_file()]
            before={p:sha(p) for p in inputs}
            proc=self.invoke(root,sample)
            self.assertEqual(proc.returncode,0,proc.stdout+proc.stderr)
            out=root/'results'/(sample+'_junction_refine_v1_1')
            s=json.loads((out/'refinement_summary.json').read_text())
            self.assertEqual(s['stats']['host_supported_pairs'],5)
            self.assertEqual(s['stats']['template_supported'],3)
            self.assertEqual(s['stats']['rescued_to_existing_template_pairs'],1)
            self.assertEqual(s['stats']['no_nearby_template'],1)
            self.assertEqual(s['stats']['no_qualifying_template'],1)
            self.assertEqual(s['stats']['supported_pairs_by_reference_compartment']['mitochondrial_reference'],1)
            self.assertTrue(all(s['validation'].values()))
            by={r['read_id']:r for r in read_table(out/'pair_refinement.tsv.gz')}
            self.assertEqual(by['pairA']['supporting_mates'],'R1;R2')
            self.assertEqual(by['pairB']['change_type'],'rescued_to_existing_template')
            self.assertEqual(by['pairC']['status'],'no_qualifying_template')
            self.assertEqual(by['pairD']['status'],'no_nearby_template')
            self.assertEqual(sum(int(r['supported_pairs']) for r in read_table(out/'identical_sequence_groups.tsv')),3)
            self.assertEqual(before,{p:sha(p) for p in inputs})

    def test_changed_pilot_table_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);sample,evidence,pilot,paths=make_runner_fixture(root)
            with (pilot/'junctions.tsv').open('a') as f:f.write('changed\n')
            proc=self.invoke(root,sample)
            self.assertNotEqual(proc.returncode,0)
            self.assertIn('hash mismatch',proc.stderr)

    def test_changed_raw_input_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);sample,evidence,pilot,paths=make_runner_fixture(root)
            with paths['R1'].open('ab') as f:f.write(b'changed')
            proc=self.invoke(root,sample)
            self.assertNotEqual(proc.returncode,0)
            self.assertIn('FASTQ changed',proc.stderr)
DONOR='TAGAGATGTTAAGCTGTTGCTTACCTTCTGAATCACCGCAACATCGGTTGAGCAGAATAGCTCCGGTCAGTGTGTCGGATCGGCGTGCCCTGAATTAGCATTAGGGCTAGTAACCCTGTACAATATGGGCCGCTCGTCCGCGTACCGGAAACCGTCTTGAACGACCAGAGACATAGTTCTAAGGTGGCTCAAGAAGGCGAGTTAAATGCACCCGCTTAGCCCTCAGAATCTGCCAATATCAATGCGACGCTGATACGCGGCTAAAGGATAGCTGTTATACAGTGCGTGCACACTCCCACGTAGGACTCCGCCCTCTTTGATGTGACCTAAAACCAAGTTGTTTCTGCTACTTAGCTGAGGGAGTTCAACCACGTGTATCATGAGTAGTAGAGCCCTCTGGATAACGGTCATCTATTGCGGGTACTAACACCCGAGTACCTTCTTGAAGTATGCTAGGCGCCCCGAAGTGTCCGAATACCCCGATTTCAACATAGCGGTGAGCGGAAAGCATCGTGATTCAGACCCAACATGCTGAGCTCAGTGATGGGTCTCCCCGTAGGTATCTATGCTGACGTCTGATCGCTCTTAGGTTTAAGGAATTTTGCAGGCTGGCTCCTGTTTATAACGGCCTGCTCTCGGAATGACAACCCAAGATACTATCGCACACGTCGTTAAGGAAGAGTAACGATGATCCTTAGGAGCTAAGTGTCAATATAGAGCTACTCTGGGCTCGATCATCCAGTTTACTGACGCGCGTTGCCCGCCGAGGACAGATTCACAGTAATGTGACTCTGTCGAAT'
