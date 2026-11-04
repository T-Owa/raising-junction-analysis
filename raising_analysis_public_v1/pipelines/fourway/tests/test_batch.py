#!/usr/bin/env python3
"""Portable batch/5-prime checks; external aligners are not executed here."""
import csv
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import textwrap
import unittest
from unittest import mock

PACKAGE = Path(__file__).resolve().parents[1]


def engine_python(source, engine='engine_v1_1', extra_args=()):
    result = subprocess.run([sys.executable, '-c', textwrap.dedent(source),
                             str(PACKAGE/engine), *map(str, extra_args)],
                            capture_output=True, text=True)
    if result.returncode:
        raise AssertionError(result.stdout + '\n' + result.stderr)
    return result.stdout


def load_batch():
    path = PACKAGE/'raising_batch.py'
    spec = importlib.util.spec_from_file_location('raising_batch_test_target', path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


# Actual runner/report code is used; only external alignment calls are mocked.
# The resulting manifest/hashes/CSVs are suitable for validator tests.
FAKE_CHILD_SCRIPT = r'''
import sys,json
from pathlib import Path
from unittest import mock
command=json.loads(sys.argv[1])
index=next(i for i,x in enumerate(command) if Path(x).name=='raising_fullrun.py')
engine=Path(command[index]).parent
sys.path.insert(0,str(engine))
import raising_fullrun as runner
import raising_threeway_classify as core
import raising_concatemer_screen as ds
args=runner.build_parser().parse_args(command[index+1:])
donor=ds.load_donor(engine/'RI_donor_800bp.fasta')
info=runner.inspect_reference(args.host_fasta,donor,args.root)
if not info['complete']:
    fa=Path(info['fasta']);fa.parent.mkdir(parents=True,exist_ok=True)
    fa.write_text('>1\nACGTACGT\n>'+core.DONOR_CONTIG+'\n'+donor+'\n')
    for suffix in ['.amb','.ann','.bwt','.pac','.sa']:
        Path(str(fa)+suffix).write_bytes(b'fixture')
    (fa.parent/'index_complete.json').write_text(json.dumps({'signature':info['signature']}))
def fake_donor(records,*args):
    for group in records.values():
        for rec in group.values():
            rec.update(dd=False,donor_mask=set(),donor_class=ds.CLASSES[2],
                       donor_details={'junction_type':'','candidate_alignments':0})
def fake_host(records,*args):
    for group in records.values():
        for rec in group.values(): rec.update(host=False,host_hits=[])
with mock.patch.object(runner.shutil,'which',return_value='/mock/tool'), \
     mock.patch.object(core,'tool_versions',return_value={'fixture':'mocked aligners'}), \
     mock.patch.object(core,'align_donor',side_effect=fake_donor), \
     mock.patch.object(core,'map_host',side_effect=fake_host):
    result=runner.run_analysis(args)
assert result['status']=='complete'
'''


def run_fixture_child(command, logpath):
    completed=subprocess.run([sys.executable,'-c',FAKE_CHILD_SCRIPT,json.dumps(list(map(str,command)))],
                             capture_output=True,text=True)
    Path(logpath).parent.mkdir(parents=True,exist_ok=True)
    Path(logpath).write_text(completed.stdout+'\n'+completed.stderr)
    return completed.returncode


class FivePrimeSequenceTests(unittest.TestCase):
    def test_full_five_prime_engine_export_with_mocked_alignments(self):
        engine_python('''
            import sys,tempfile,json,gzip
            from pathlib import Path
            from unittest import mock
            sys.path.insert(0,sys.argv[1])
            import raising_fullrun as runner
            import raising_fullrun_report as report
            import raising_threeway_classify as c
            import raising_concatemer_screen as d
            donor=d.load_donor(Path(sys.argv[1])/'RI_donor_800bp.fasta')
            rc=d.revcomp(donor); start=rc.index(d.P3['5p'])
            short=rc[start:start+55]+'A'*20+c.R1_ADAPTER
            ddseq=d.P3['5p']+d.revcomp(donor[:65])+d.revcomp(donor[750:])
            with tempfile.TemporaryDirectory() as td:
                root=Path(td); host=root/'host.fa'; host.write_text('>1\\nACGTACGT\\n')
                raw=root/'vendor_results/Raw_fastq'; raw.mkdir(parents=True)
                sample='fixture_5prm_S1'
                for mate in ['R1','R2']:
                    with gzip.open(raw/f'{sample}_L001_{mate}_001.fastq.gz','wt') as f:
                        for rid in ['short','host','dd']:
                            seq=(ddseq if rid=='dd' else short) if mate=='R1' else 'ACGT'*30
                            f.write(f'@{rid}/{mate[-1]}\\n{seq}\\n+\\n'+('I'*len(seq))+'\\n')
                info=runner.inspect_reference(host,donor,root)
                fa=Path(info['fasta']); fa.parent.mkdir(parents=True)
                fa.write_text('>1\\nACGTACGT\\n>'+c.DONOR_CONTIG+'\\n'+donor+'\\n')
                for suffix in ['.amb','.ann','.bwt','.pac','.sa']:
                    Path(str(fa)+suffix).write_bytes(b'fixture')
                (fa.parent/'index_complete.json').write_text(json.dumps({'signature':info['signature']}))
                def fake_donor(records,donor,sample,threads,work,commands):
                    assert '_5prm' in sample
                    for mate,group in records.items():
                        for rid,rec in group.items():
                            if mate=='R1': assert rec['raw_start0']==23
                            rec.update(dd=(rid=='dd' and mate=='R1'),donor_mask=set(),
                                donor_class=d.CLASSES[0] if rid=='dd' else d.CLASSES[2],
                                donor_details={'unaligned_gap_bp':0,'alignment_overlap_bp':0,
                                  'junction_type':'H-T_like','candidate_alignments':1} if rid=='dd' and mate=='R1'
                                  else {'junction_type':'','candidate_alignments':0})
                def fake_host(records,*args):
                    for mate,group in records.items():
                        for rid,rec in group.items():
                            rec['host']=rid=='host' and mate=='R2'; rec['host_hits']=[]
                argv=['--root',str(root),'--sample',sample,'--host-fasta',str(host),'--chunk-pairs','2']
                with mock.patch.object(runner.shutil,'which',return_value='/mock/tool'), \\
                     mock.patch.object(c,'tool_versions',return_value={'fixture':'mocked aligners'}), \\
                     mock.patch.object(c,'align_donor',side_effect=fake_donor), \\
                     mock.patch.object(c,'map_host',side_effect=fake_host):
                    result=runner.run_analysis(runner.build_parser().parse_args(argv))
                assert result['status']=='complete'
                assert result['side']=='5p' and result['config']['side']=='5p'
                assert result['analysis_version']=='1.1.0'
                assert result['class_counts']==dict(zip(report.CLASSES,[1,1,1,0]))
                assert all(result['validation'].values())
                out=root/'results'/f'{sample}_fourway_full_v1_1'
                methods=(out/'methods.txt').read_text()
                assert 'Library direction: 5p' in methods
                assert '3-prime libraries only' not in methods
        ''')

    def test_five_prime_reference_orientation_and_mask(self):
        engine_python('''
            import sys
            from pathlib import Path
            sys.path.insert(0,sys.argv[1])
            import raising_concatemer_screen as d
            donor=d.load_donor(Path(sys.argv[1])/'RI_donor_800bp.fasta')
            oriented=d.revcomp(donor)
            assert oriented.index(d.P3['5p'])==712
            assert donor[65:88]==d.revcomp(d.P3['5p'])
            mask=d.primer_reference_mask(donor,'fixture_5prm_S1')
            assert mask==set(range(65,88))
            seq=d.revcomp(donor[48:88])
            fields=['U1','1','40','88','49','100','40','40','0','0','1e-20','80',seq,seq]
            h=d.parse_hit(fields,{'U1':seq},donor,mask)
            assert h.strand=='-'
            assert len(h.matched_query)==40
            assert len(h.independent_query)==17, len(h.independent_query)
        ''')

    def test_five_prime_outward_short_pattern_and_adapter(self):
        engine_python('''
            import sys
            from pathlib import Path
            sys.path.insert(0,sys.argv[1])
            import raising_fullrun_report as r
            import raising_threeway_classify as c
            import raising_concatemer_screen as d
            donor=d.load_donor(Path(sys.argv[1])/'RI_donor_800bp.fasta')
            rc=d.revcomp(donor); start=rc.index(d.P3['5p'])
            prefix=rc[start:start+55]
            assert prefix[-1]!='A'
            seq=prefix+'A'*20+c.R1_ADAPTER
            pat=r.short_product_pattern(seq,donor,'5p')
            assert pat['short_pattern'] and pat['adapter_confirmed'] and pat['adapter_exact']
            assert pat['prefix_bp']==55
            assert pat['bases_beyond_P3']==32
            assert pat['prefix_donor_oriented_start_1based']==713
            assert not r.short_product_pattern(seq,donor,'3p')['short_pattern']
            rec=c.clean_record(('@x',seq,'+','I'*len(seq)),'R1','5p')
            assert rec['raw_start0']==23
            assert rec['sequence']==prefix[23:]
            assert 'confirmed_polyA_adapter_trimmed' in rec['flags']
        ''')

    def test_five_prime_primer_only_and_unconfirmed_polya_are_distinct(self):
        engine_python('''
            import sys
            from pathlib import Path
            sys.path.insert(0,sys.argv[1])
            import raising_fullrun_report as r
            import raising_threeway_classify as c
            import raising_concatemer_screen as d
            donor=d.load_donor(Path(sys.argv[1])/'RI_donor_800bp.fasta')
            p=d.P3['5p']
            raw=p+'A'*20+c.R1_ADAPTER
            result=r.short_product_pattern(raw,donor,'5p')
            assert result['short_pattern'] and result['primer_only_prefix']
            assert result['bases_beyond_P3']==0
            raw=p+'A'*20+'CCGTCGTCGT'
            result=r.short_product_pattern(raw,donor,'5p')
            assert result['short_pattern'] and not result['adapter_confirmed']
            rec=c.clean_record(('@x',raw,'+','I'*len(raw)),'R1','5p')
            assert rec['technical_suffix_bp']==0
            assert rec['sequence'].startswith('A'*20)
        ''')

    def test_reverse_strand_concatemer_anchors_remain_supported(self):
        engine_python('''
            import sys
            from pathlib import Path
            from types import SimpleNamespace
            sys.path.insert(0,sys.argv[1])
            import raising_concatemer_screen as d
            donor=d.load_donor(Path(sys.argv[1])/'RI_donor_800bp.fasta')
            a=d.revcomp(donor[:65]); b=d.revcomp(donor[750:800]); seq=a+b
            mask=d.primer_reference_mask(donor,'fixture_5prm_S1')
            f1=['U1','1','65','65','1','100','65','65','0','0','1e-20','130',a,a]
            f2=['U1','66','115','800','751','100','50','50','0','0','1e-20','100',b,b]
            hs=[d.parse_hit(f,{'U1':seq},donor,mask) for f in (f1,f2)]
            p=SimpleNamespace(min_anchor=30,min_identity=95,max_gap=50,max_overlap=20,min_rewind=20)
            result=d.classify_hits(hs,len(seq),p)
            assert result['classification']==d.CLASSES[0]
            assert result['best']['kind']=='H-T_like'
            assert result['best']['left_independent_matches']==65
            assert result['best']['right_independent_matches']==50
        ''')

    def test_r2_support_overrides_five_prime_short_r1(self):
        engine_python('''
            import sys
            from pathlib import Path
            sys.path.insert(0,sys.argv[1])
            import raising_fullrun_report as r
            import raising_threeway_classify as c
            import raising_concatemer_screen as d
            donor=d.load_donor(Path(sys.argv[1])/'RI_donor_800bp.fasta')
            raw=d.P3['5p']+'A'*20+c.R1_ADAPTER
            pat=r.short_product_pattern(raw,donor,'5p')
            def rec(dd=False,host=False):
                return dict(dd=dd,host=host,sequence='ACGT'*30,donor_mask=set(),flags=set(),host_hits=[],
                    donor_details={'unaligned_gap_bp':0,'junction_type':'H-T_like'} if dd else {})
            assert r.classify_pair([rec(),rec(host=True)],pat)['classification']==r.CLASSES[1]
            both=r.classify_pair([rec(),rec(dd=True,host=True)],pat)
            assert both['classification']==r.CLASSES[0] and both['both_supported']==1
        ''')

    def test_engine_side_detection_accepts_5p_and_rejects_ambiguity(self):
        engine_python('''
            import sys
            sys.path.insert(0,sys.argv[1])
            import raising_fullrun as r
            # Side validation precedes missing-file checks. A valid 5p sample
            # must not fail with the old 3-prime-only error.
            args=r.build_parser().parse_args(['--sample','fixture_5prm_S1','--host-fasta','/missing/host.fa'])
            try: r.run_analysis(args)
            except ValueError as e: assert 'Host FASTA not found' in str(e),str(e)
            else: raise AssertionError('Missing host should fail')
            for name in ['fixture_S1','fixture_3prm_5prm_S1']:
                args=r.build_parser().parse_args(['--sample',name,'--host-fasta','/missing/host.fa'])
                try: r.run_analysis(args)
                except ValueError as e: assert 'Host FASTA not found' not in str(e),str(e)
                else: raise AssertionError('Invalid side should fail')
        ''')


class DiscoveryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.batch=load_batch()

    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.path=Path(self.temp.name)

    def tearDown(self): self.temp.cleanup()

    def put(self,sample,lane='L001',mate='R1',part='001'):
        path=self.path/f'{sample}_{lane}_{mate}_{part}.fastq.gz'
        path.write_bytes(b'fixture')
        return path

    def pair(self,sample,lane='L001',part='001'):
        return [self.put(sample,lane,m,part) for m in ('R1','R2')]

    def test_both_directions_and_multilane_are_grouped_once(self):
        self.pair('sample001_3prm_S1','L002'); self.pair('sample001_3prm_S1','L001')
        self.pair('sample001_5prm_S2','L001')
        libs=self.batch.discover_libraries(self.path)
        self.assertEqual(len(libs),2)
        by={x['sample']:x for x in libs}
        self.assertEqual(by['sample001_3prm_S1']['side'],'3p')
        self.assertEqual(by['sample001_5prm_S2']['side'],'5p')
        self.assertEqual(len(by['sample001_3prm_S1']['r1']),2)
        self.assertTrue(str(by['sample001_3prm_S1']['r1'][0]).endswith('_L001_R1_001.fastq.gz'))

    def test_missing_r2_stops_discovery(self):
        self.put('sample001_3prm_S1')
        with self.assertRaises(ValueError): self.batch.discover_libraries(self.path)

    def test_orphan_r2_stops_discovery(self):
        self.put('sample001_3prm_S1',mate='R2')
        with self.assertRaises(ValueError): self.batch.discover_libraries(self.path)

    def test_unknown_or_ambiguous_side_stops_discovery(self):
        for name in ('sample001_S1','sample001_3prm_5prm_S1','3prm_foo',
                     'X_3prm_Y_5prmBAD','X_3prm_3prmBAD','X_5prm_3prmBAD'):
            with self.subTest(name=name):
                for old in self.path.iterdir(): old.unlink()
                self.pair(name)
                with self.assertRaises(ValueError): self.batch.discover_libraries(self.path)

    def test_empty_directory_has_no_libraries(self):
        with self.assertRaises(ValueError): self.batch.discover_libraries(self.path)


class FrozenCoreTests(unittest.TestCase):

    def test_alignment_core_is_identical_between_engines(self):
        for name in ('raising_threeway_classify.py','raising_concatemer_screen.py','RI_donor_800bp.fasta'):
            self.assertEqual((PACKAGE/'legacy_v1'/name).read_bytes(),(PACKAGE/'engine_v1_1'/name).read_bytes(),name)



class BatchIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls): cls.batch=load_batch()

    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.root=Path(self.temp.name)
        self.raw=self.root/'vendor_results/Raw_fastq'; self.raw.mkdir(parents=True)
        self.host=self.root/'host.fa'; self.host.write_text('>1\nACGTACGT\n')
        self.out=self.root/'results/all_libraries_fourway_batch_v1'
        donor=''.join(x for x in (PACKAGE/'legacy_v1/RI_donor_800bp.fasta').read_text().splitlines() if not x.startswith('>'))
        rc=donor.translate(str.maketrans('ACGT','TGCA'))[::-1]
        prefixes={'3p':donor[718:773], '5p':rc[712:767]}
        self.samples={'3p':'sample001_3prm_S1','5p':'sample001_5prm_S2'}
        for side,sample in self.samples.items():
            for mate in ('R1','R2'):
                with gzip.open(self.raw/f'{sample}_L001_{mate}_001.fastq.gz','wt') as f:
                    for i in (1,2):
                        seq=(prefixes[side]+'A'*20+'TACTGCTTGACTGACCTGCTGTC') if i==1 and mate=='R1' else 'CGTC'*30
                        f.write(f'@{sample}:{i}/{mate[-1]}\n{seq}\n+\n'+('I'*len(seq))+'\n')
        self.argv=['--root',str(self.root),'--host-fasta',str(self.host),'--chunk-pairs','1']
        self.calls=[]
        self.failure_sample=None
        self.interrupt_sample=None
        def child(command,logpath):
            sample=command[command.index('--sample')+1]
            self.calls.append(list(command))
            if sample==self.interrupt_sample: raise KeyboardInterrupt()
            if sample==self.failure_sample:
                Path(logpath).write_text('Intentional fixture failure')
                return 1
            return run_fixture_child(command,logpath)
        self.patches=[mock.patch.object(self.batch.shutil,'which',return_value='/mock/tool'),
                      mock.patch.object(self.batch,'run_child',side_effect=child)]
        self.mocks=[p.start() for p in self.patches]

    def tearDown(self):
        for patch in reversed(self.patches): patch.stop()
        self.temp.cleanup()

    def run_batch(self,extra=()):
        return self.batch.run_batch(self.batch.build_parser().parse_args(self.argv+list(extra)))

    def rows(self,name):
        with (self.out/name).open(encoding='utf-8-sig') as f: return list(csv.DictReader(f))

    def test_both_sides_full_export_and_completed_resume(self):
        result=self.run_batch()
        self.assertEqual(result['status'],'complete')
        self.assertEqual(result['completed_libraries'],2)
        self.assertEqual(len(self.rows('all_figure_counts.csv')),8)
        wide=self.rows('all_figure_wide.csv')
        self.assertEqual({r['library_side'] for r in wide},{'3p','5p'})
        self.assertEqual({r['analysis_version'] for r in wide},{'1.0.0','1.1.0'})
        self.assertEqual({int(r['total_pairs']) for r in wide},{2})
        self.assertEqual({r['mode'] for r in self.rows('library_status.csv')},{'new'})
        repeat=self.run_batch()
        self.assertEqual(repeat['status'],'complete')
        self.assertEqual({r['mode'] for r in self.rows('library_status.csv')},{'reused_completed'})
        self.assertTrue(all('--resume' in c for c in self.calls[-2:]))
        self.assertTrue(any((self.out/'previous_tables').rglob('all_figure_counts.csv')))

    def test_index_build_permission_is_forwarded_to_each_engine(self):
        result=self.run_batch(['--allow-index-build'])
        self.assertEqual(result['status'],'complete')
        self.assertEqual(len(self.calls),2)
        self.assertTrue(all('--allow-index-build' in c for c in self.calls))

    def test_failed_library_is_not_silently_omitted_from_success(self):
        self.failure_sample=self.samples['5p']
        result=self.run_batch()
        self.assertEqual(result['status'],'partial')
        self.assertEqual(result['completed_libraries'],1)
        self.assertEqual(result['failed_libraries'],1)
        self.assertFalse((self.out/'all_figure_counts.csv').exists())
        self.assertEqual(len(self.rows('partial_figure_counts.csv')),4)
        status={r['sample']:r['status'] for r in self.rows('library_status.csv')}
        self.assertEqual(status[self.samples['3p']],'complete')
        self.assertEqual(status[self.samples['5p']],'failed')
        self.failure_sample=None
        self.assertEqual(self.run_batch()['status'],'complete')
        modes={r['sample']:r['mode'] for r in self.rows('library_status.csv')}
        self.assertEqual(modes[self.samples['3p']],'reused_completed')
        self.assertEqual(modes[self.samples['5p']],'new')

    def test_interrupt_preserves_completed_and_reports_pending(self):
        self.interrupt_sample=self.samples['5p']
        result=self.run_batch()
        self.assertEqual(result['status'],'interrupted')
        self.assertEqual(result['completed_libraries'],1)
        self.assertFalse((self.out/'all_figure_counts.csv').exists())
        saved=json.loads((self.out/'batch_status.json').read_text())
        self.assertEqual(saved['status'],'interrupted')
        self.interrupt_sample=None
        self.assertEqual(self.run_batch()['status'],'complete')

    def test_plan_only_does_not_launch_any_child(self):
        result=self.run_batch(['--plan-only'])
        self.assertEqual(result['status'],'plan_only')
        self.assertEqual(self.calls,[])
        self.assertEqual(len(self.rows('library_plan.csv')),2)
        self.assertFalse((self.out/'all_figure_counts.csv').exists())
        self.assertEqual(self.run_batch()['status'],'complete')

    def test_selection_change_requires_separate_batch_output(self):
        self.run_batch(['--plan-only'])
        with self.assertRaisesRegex(ValueError,'selection|selection/reference'):
            self.run_batch(['--samples',self.samples['3p']])

    def test_corrupted_completed_evidence_is_rejected_on_resume(self):
        self.run_batch()
        evidence=self.root/'results'/f"{self.samples['3p']}_fourway_full_v1"/'chunks/chunk_000001/pair_classification.tsv'
        evidence.write_text(evidence.read_text()+'\n')
        result=self.run_batch()
        self.assertEqual(result['status'],'partial')
        self.assertEqual(result['failed_libraries'],1)
        self.assertFalse((self.out/'all_figure_counts.csv').exists())

    def test_validator_detects_tampered_figure_csv(self):
        self.run_batch()
        lib=self.batch.discover_libraries(self.raw)[0]
        resultdir=self.root/'results'/f"{lib['sample']}_fourway_full_v1"
        path=resultdir/'figure_counts.csv'
        path.write_text(path.read_text(encoding='utf-8-sig').replace(',0,',',1,',1),encoding='utf-8-sig')
        with self.assertRaisesRegex(ValueError,'checksum'):
            self.batch.validate_result(lib['sample'],lib['side'],resultdir,{},lib)

    def test_validator_checks_count_arithmetic_even_if_checksum_updated(self):
        self.run_batch()
        lib=self.batch.discover_libraries(self.raw)[0]
        resultdir=self.root/'results'/f"{lib['sample']}_fourway_full_v1"
        path=resultdir/'figure_counts.csv'
        with path.open(encoding='utf-8-sig') as f:
            reader=csv.DictReader(f); fields=reader.fieldnames; rows=list(reader)
        rows[0]['percent']='99'
        with path.open('w',encoding='utf-8-sig',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(rows)
        summarypath=resultdir/'run_summary.json'; summary=json.loads(summarypath.read_text())
        summary['final_output_sha256']['figure_counts.csv']=hashlib.sha256(path.read_bytes()).hexdigest()
        summarypath.write_text(json.dumps(summary))
        with self.assertRaisesRegex(ValueError,'percentage'):
            self.batch.validate_result(lib['sample'],lib['side'],resultdir,{},lib)


if __name__=='__main__': unittest.main()
