#!/usr/bin/env python3
"""Sequential, resumable RAISING full-input analysis of all 3'/5' libraries."""
import argparse
import collections
import contextlib
import csv
import datetime
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
import uuid

VERSION = '1.0.0'
BASE = Path(__file__).resolve().parent
CLASSES = ('donor_rearrangement_supported','host_mappable_other',
           'short_donor_primer_candidate','other_unresolved')
SOURCE_FILES = ('raising_fullrun.py','raising_fullrun_report.py',
                'raising_threeway_classify.py','raising_concatemer_screen.py')
THRESHOLDS = {'donor_anchor_independent_matches':30,'donor_min_identity_percent':95,
    'donor_max_gap_bp':50,'donor_max_overlap_bp':20,'donor_min_rewind_bp':20,
    'host_independent_matches':30,'host_min_identity_percent':95,
    'host_min_mapq':20,'host_min_base_quality':20,
    'donor_base_quality_filter':'none_in_preserved_core'}
FASTQ_PATTERN = re.compile(r'^(?P<sample>.+?)(?P<lane>_L\d{3})?_R(?P<mate>[12])_(?P<part>\d{3})\.fastq(?:\.gz)?$')


def parse_side(sample):
    if not re.fullmatch(r'[A-Za-z0-9._-]+',sample) or sample in ('.','..'):
        raise ValueError('Invalid sample name: '+sample)
    tokens = re.findall(r'_([35])prm(?=_|$)',sample)
    if (len(tokens)!=1 or ('_3prm' in sample and '_5prm' in sample)
            or re.search(r'_[35]prm(?!_|$)',sample)):
        raise ValueError('Sample must contain exactly one _3prm or _5prm token: '+sample)
    return tokens[0]+'p'


def discover_libraries(rawdir):
    rawdir = Path(rawdir).expanduser().resolve()
    if not rawdir.is_dir():
        raise ValueError('Raw FASTQ directory not found: '+str(rawdir))
    groups = collections.defaultdict(dict)
    unknown = []
    for path in sorted(rawdir.iterdir()):
        if not path.is_file() or not (path.name.endswith('.fastq.gz') or path.name.endswith('.fastq')):
            continue
        m = FASTQ_PATTERN.fullmatch(path.name)
        if not m:
            unknown.append(path.name); continue
        sample = m['sample']; parse_side(sample)
        key = (m['lane'] or '',m['part'],m['mate'])
        if key in groups[sample]:
            raise ValueError('Duplicate lane/part/mate, including compressed and uncompressed copies: '+sample)
        groups[sample][key] = path.resolve()
    if unknown:
        raise ValueError('Unrecognized FASTQ names (none silently omitted): '+', '.join(unknown))
    if not groups:
        raise ValueError('No paired RAISING FASTQ libraries were found')
    libraries = []
    for sample, files in sorted(groups.items()):
        partitions = sorted({(lane,part) for lane,part,mate in files})
        r1, r2 = [], []
        for lane,part in partitions:
            if (lane,part,'1') not in files or (lane,part,'2') not in files:
                raise ValueError(f'Missing R1/R2 partner: {sample}{lane}, part {part}')
            r1.append(files[lane,part,'1']); r2.append(files[lane,part,'2'])
        libraries.append({'sample':sample,'side':parse_side(sample),'r1':r1,'r2':r2})
    return libraries


def engine_info(side):
    if side=='3p':
        return BASE/'legacy_v1'/'raising_fullrun.py','1.0.0','_fourway_full_v1'
    if side=='5p':
        return BASE/'engine_v1_1'/'raising_fullrun.py','1.1.0','_fourway_full_v1_1'
    raise ValueError('Unknown library side')


def sha256_file(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''): h.update(block)
    return h.hexdigest()


def atomic_json(path,data):
    path = Path(path); tmp = path.with_name(path.name+'.tmp')
    with tmp.open('w',encoding='utf-8') as f:
        json.dump(data,f,ensure_ascii=False,indent=2,default=str)
        f.write('\n'); f.flush(); os.fsync(f.fileno())
    tmp.replace(path)


def read_csv(path):
    with Path(path).open(encoding='utf-8-sig',newline='') as f:
        return list(csv.DictReader(f))


def write_csv(path,rows,fields=None):
    path = Path(path); tmp = path.with_name(path.name+'.tmp')
    if fields is None:
        fields = list(dict.fromkeys(k for row in rows for k in row))
    with tmp.open('w',encoding='utf-8-sig',newline='') as f:
        w = csv.DictWriter(f,fieldnames=fields)
        w.writeheader(); w.writerows(rows)
    tmp.replace(path)


def _integer(value,label):
    if isinstance(value,bool): raise ValueError('Boolean count: '+label)
    n = int(value)
    if str(n)!=str(value) or n<0: raise ValueError('Invalid nonnegative integer: '+label)
    return n


def validate_result(sample,side,outdir,expected_config,expected_inputs):
    """Validate a child-completed result and its source tables before merging."""
    outdir = Path(outdir)
    summary = json.loads((outdir/'run_summary.json').read_text())
    if summary.get('status')!='complete' or summary.get('classification_available') is not True:
        raise ValueError('Analysis is not complete: '+sample)
    if summary.get('sample')!=sample or summary.get('side')!=side:
        raise ValueError('Sample/direction mismatch: '+sample)
    ref = summary.get('reference',{})
    for digest in (ref.get('host_fasta_sha256',''),ref.get('signature',{}).get('donor_sha256','')):
        if not re.fullmatch('[0-9a-f]{64}',digest):
            raise ValueError('Reference identity checksum missing: '+sample)
    if not summary.get('tool_versions'):
        raise ValueError('Software version provenance missing: '+sample)
    validation = summary.get('validation',{})
    if not validation or not all(v is True for v in validation.values()):
        raise ValueError('Analysis validation did not pass: '+sample)
    for key,value in expected_config.items():
        observed = summary.get('source_sha256',{}) if key=='_source_sha256' else summary.get('config',{}).get(key)
        if observed!=value:
            raise ValueError('Analysis configuration/source mismatch: '+sample+' '+key)
    lanes = summary.get('inputs',{}).get('lanes',[])
    if len(lanes)!=len(expected_inputs['r1']) or len(lanes)!=len(expected_inputs['r2']):
        raise ValueError('Input lane count mismatch: '+sample)
    for i,lane in enumerate(lanes):
        for mate in ('R1','R2'):
            path = Path(expected_inputs[mate.lower()][i]).resolve()
            rec = lane.get(mate,{})
            st = path.stat()
            if rec.get('path')!=str(path) or rec.get('size')!=st.st_size or rec.get('mtime_ns')!=st.st_mtime_ns:
                raise ValueError('Input file changed or differs from analyzed file: '+str(path))
            if not re.fullmatch('[0-9a-f]{64}',rec.get('sha256','')):
                raise ValueError('Input SHA-256 missing: '+sample)
    n = _integer(summary['total_pairs'],'total_pairs')
    counts = summary['class_counts']
    if n==0 or set(counts)!=set(CLASSES) or sum(_integer(counts[c],c) for c in CLASSES)!=n:
        raise ValueError('Invalid exclusive class counts: '+sample)
    files = {'figure_counts':'figure_counts.csv','figure_wide':'figure_wide.csv',
             'supplementary_counts':'supplementary_counts.csv'}
    rows = {}
    for key,name in files.items():
        path = outdir/name
        if sha256_file(path)!=summary.get('final_output_sha256',{}).get(name):
            raise ValueError('Figure/source table checksum mismatch: '+str(path))
        rows[key] = read_csv(path)
        for r in rows[key]:
            if r.get('sample')!=sample or r.get('library_side')!=side:
                raise ValueError('CSV sample/direction mismatch: '+str(path))
    primary = rows['figure_counts']
    if len(primary)!=4 or {r['category'] for r in primary}!=set(CLASSES):
        raise ValueError('Expected exactly four exclusive class rows: '+sample)
    for r in primary:
        count = _integer(r['count'],'count')
        if count!=counts[r['category']] or _integer(r['total_pairs'],'total_pairs')!=n:
            raise ValueError('CSV count differs from completed summary: '+sample)
        if r.get('denominator')!='all_input_read_pairs':
            raise ValueError('Unexpected denominator: '+sample)
        if not math.isclose(float(r['fraction']),count/n,rel_tol=0,abs_tol=1e-12) or not math.isclose(float(r['percent']),100*count/n,rel_tol=0,abs_tol=1e-9):
            raise ValueError('CSV percentage differs from counts: '+sample)
    wide = rows['figure_wide']
    if len(wide)!=1 or _integer(wide[0]['total_pairs'],'total_pairs')!=n:
        raise ValueError('Expected one wide row: '+sample)
    for c in CLASSES:
        if _integer(wide[0][c+'_count'],c)!=counts[c] or not math.isclose(float(wide[0][c+'_percent']),100*counts[c]/n,rel_tol=0,abs_tol=1e-9):
            raise ValueError('Wide table differs from summary: '+sample)
    for r in rows['supplementary_counts']:
        count = _integer(r['count'],'supplementary count')
        den = _integer(r['denominator_count'],'supplementary denominator')
        expected_den = n if r['scope']=='all_input_pairs_nonexclusive' else counts[CLASSES[2]] if r['scope']=='short_donor_primer_candidate_class' else None
        if den!=expected_den or count>den:
            raise ValueError('Invalid supplementary denominator/count: '+sample)
        if den and not math.isclose(float(r['percent']),100*count/den,rel_tol=0,abs_tol=1e-9):
            raise ValueError('Invalid supplementary percentage: '+sample)
        if not den and r['percent']!='':
            raise ValueError('Zero-denominator percentage must be blank: '+sample)
    return {'summary':summary,**rows,'output_dir':str(outdir)}


def aggregate_results(libraries,statuses,outdir):
    outdir = Path(outdir); outdir.mkdir(parents=True,exist_ok=True)
    expected = {x['sample']:x for x in libraries}
    if len(expected)!=len(libraries): raise ValueError('Duplicate planned sample')
    by_id = {x['sample']:x for x in statuses}
    if len(by_id)!=len(statuses) or set(by_id)!=set(expected):
        raise ValueError('Status rows do not match planned libraries')
    success = [by_id[x['sample']] for x in libraries if by_id[x['sample']]['status']=='complete']
    all_complete = len(success)==len(libraries)
    status_rows = []
    for lib in libraries:
        s = by_id[lib['sample']]
        meta = s.get('validated',{}).get('summary',{})
        status_rows.append({'sample':lib['sample'],'library_side':lib['side'],
            'status':s['status'],'mode':s.get('mode',''),'total_pairs':meta.get('total_pairs',''),
            'analysis_version':meta.get('analysis_version',''),
            'elapsed_seconds':s.get('elapsed_seconds',''),'result_directory':s.get('output_dir',''),
            'log_file':s.get('log_file',''),'error':s.get('error','')})
    write_csv(outdir/'library_status.csv',status_rows)
    # Cross-library results must describe the same reference and alignment criteria.
    comparability = set()
    for s in success:
        v = s['validated']; meta = v['summary']
        if meta['sample']!=s['sample'] or meta['side']!=expected[s['sample']]['side']:
            raise ValueError('Validated result attached to wrong library')
        comparability.add(json.dumps({'host_sha256':meta.get('reference',{}).get('host_fasta_sha256'),
            'donor_sha256':meta.get('reference',{}).get('signature',{}).get('donor_sha256'),
            'thresholds':{k:meta['config'].get(k) for k in THRESHOLDS},
            'tools':meta.get('tool_versions',{})},sort_keys=True))
    if len(comparability)>1:
        raise ValueError('Completed libraries use different reference, thresholds, or tool versions; do not merge')
    prefix = 'all' if all_complete else 'partial'
    outputs = []
    for key in ('figure_counts','figure_wide','supplementary_counts'):
        combined = []
        for s in success:
            for row in s['validated'][key]:
                combined.append({**row,'analysis_version':s['validated']['summary']['analysis_version'],
                                 'result_directory':str(s['output_dir'])})
        if combined:
            path = outdir/(prefix+'_'+key+'.csv')
            write_csv(path,combined); outputs.append(path.name)
    return {'status':'complete' if all_complete else 'partial',
            'total_libraries':len(libraries),'completed_libraries':len(success),
            'failed_libraries':sum(s['status']=='failed' for s in statuses),
            'pending_libraries':sum(s['status'] in ('pending','running','interrupted') for s in statuses),
            'aggregate_files':outputs+['library_status.csv']}


def run_child(command,logpath):
    with Path(logpath).open('w',encoding='utf-8') as log:
        log.write('COMMAND '+json.dumps(command)+'\n'); log.flush()
        proc = subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                                text=True,encoding='utf-8',errors='replace')
        try:
            for line in proc.stdout:
                print(line,end='',flush=True); log.write(line); log.flush()
            return proc.wait()
        except BaseException:
            if proc.poll() is None:
                proc.terminate()
                try: proc.wait(timeout=5)
                except subprocess.TimeoutExpired: proc.kill(); proc.wait()
            raise


@contextlib.contextmanager
def lock_output(outdir):
    with (outdir/'.batch.lock').open('a+') as f:
        try: fcntl.flock(f.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError as e: raise ValueError('Another batch is using this output directory') from e
        try: yield
        finally: fcntl.flock(f.fileno(),fcntl.LOCK_UN)


def checkpoint(outdir,manifest,statuses,state,started,extra=None):
    data = {**manifest,'status':state,'elapsed_seconds':time.perf_counter()-started,
            'libraries':[{k:v for k,v in s.items() if k!='validated'} for s in statuses]}
    if extra: data.update(extra)
    atomic_json(outdir/'batch_status.json',data)


def build_parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path.home()/'raising')
    p.add_argument('--host-fasta',type=Path,required=True)
    p.add_argument('--threads',type=int,default=4)
    p.add_argument('--chunk-pairs',type=int,default=10000)
    p.add_argument('--outdir',type=Path)
    p.add_argument('--samples',nargs='+')
    p.add_argument('--plan-only',action='store_true')
    p.add_argument('--allow-index-build',action='store_true',
                   help='Allow the existing engines to build a missing human+donor BWA index (needed for a fresh root)')
    return p


def run_batch(args):
    started = time.perf_counter()
    root = args.root.expanduser().resolve(); host = args.host_fasta.expanduser().resolve()
    if not host.is_file(): raise ValueError('Host FASTA not found: '+str(host))
    if args.threads<1 or args.chunk_pairs<1: raise ValueError('Threads/chunk size must be positive')
    libraries = discover_libraries(root/'vendor_results'/'Raw_fastq')
    if args.samples:
        wanted = set(args.samples)
        if len(wanted)!=len(args.samples): raise ValueError('Repeated --samples entry')
        missing = wanted-{x['sample'] for x in libraries}
        if missing: raise ValueError('Requested samples not found: '+', '.join(sorted(missing)))
        libraries = [x for x in libraries if x['sample'] in wanted]
    outdir = (args.outdir or root/'results'/'all_libraries_fourway_batch_v1').expanduser().resolve()
    for path in [host]+[p for x in libraries for p in x['r1']+x['r2']]:
        if path==outdir or outdir in path.parents: raise ValueError('Input file is inside batch output directory')
    counts = collections.Counter(x['side'] for x in libraries)
    print(f"Found {len(libraries)} libraries: 3p={counts['3p']}, 5p={counts['5p']}. Sequential full-input processing.",flush=True)
    for i,lib in enumerate(libraries,1): print(f"  {i:02d} {lib['sample']} ({lib['side']}, {len(lib['r1'])} lane/part pair(s))")
    if not args.plan_only:
        for tool in ('bwa','blastn','makeblastdb'):
            if not shutil.which(tool): raise ValueError(tool+' not found; activate the raising environment')
    outdir.mkdir(parents=True,exist_ok=True)
    with lock_output(outdir):
        manifest_path = outdir/'batch_manifest.json'
        signature = {'batch_version':VERSION,'root':str(root),'host_fasta':str(host),
            'libraries':[{'sample':x['sample'],'side':x['side'],'r1':list(map(str,x['r1'])),'r2':list(map(str,x['r2']))} for x in libraries]}
        if manifest_path.exists():
            old = json.loads(manifest_path.read_text())
            if old['selection']!=signature:
                raise ValueError('Batch selection/reference changed. Use a different --outdir to keep results separate.')
        else:
            existing = {p.name for p in outdir.iterdir()}-{'.batch.lock'}
            if existing: raise ValueError('Unrecognized files in batch output directory: '+', '.join(existing))
            atomic_json(manifest_path,{'selection':signature})
        write_csv(outdir/'library_plan.csv',[{'sample':x['sample'],'library_side':x['side'],
            'R1_files':';'.join(map(str,x['r1'])),'R2_files':';'.join(map(str,x['r2']))} for x in libraries])
        if args.plan_only:
            print('Plan saved: '+str(outdir/'library_plan.csv')+'; no mapping was run.')
            return {'status':'plan_only','total_libraries':len(libraries)}
        run_id = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'_'+uuid.uuid4().hex[:8]
        # Preserve previous batch tables; avoid presenting stale final tables as a new run.
        previous = [p for p in outdir.iterdir() if p.name in ('library_status.csv',) or
                    (p.name.startswith(('all_','partial_')) and p.suffix=='.csv')]
        if previous:
            archive = outdir/'previous_tables'/run_id; archive.mkdir(parents=True)
            for path in previous: path.rename(archive/path.name)
        logdir = outdir/'logs'/run_id; logdir.mkdir(parents=True)
        manifest = {'run_id':run_id,'batch_version':VERSION,'total_libraries':len(libraries),
                    'batch_script_sha256':sha256_file(Path(__file__)),
                    'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
        statuses = [{'sample':x['sample'],'side':x['side'],'status':'pending','error':''} for x in libraries]
        interrupted = False
        checkpoint(outdir,manifest,statuses,'running',started)
        for i,(lib,s) in enumerate(zip(libraries,statuses),1):
            begin = time.perf_counter()
            engine,version,suffix = engine_info(lib['side'])
            result_dir = root/'results'/(lib['sample']+suffix)
            s.update(output_dir=str(result_dir),log_file=str(logdir/(lib['sample']+'.log')),status='running')
            try:
                # Batch outputs must not contain or coincide with per-library result folders.
                if outdir==result_dir or outdir in result_dir.parents or result_dir in outdir.parents:
                    raise ValueError('Batch and per-library output directories must be separate')
                previous_manifest = result_dir/'run_manifest.json'
                resume = previous_manifest.is_file()
                cfg = {}
                prior_complete = False
                if resume:
                    prior = json.loads(previous_manifest.read_text())
                    cfg = prior['signature_data']['config']
                    if (result_dir/'run_summary.json').is_file():
                        prior_complete = json.loads((result_dir/'run_summary.json').read_text()).get('status')=='complete'
                s['mode'] = 'reused_completed' if prior_complete else 'resumed' if resume else 'new'
                threads = cfg.get('threads',args.threads); chunk = cfg.get('chunk_pairs',args.chunk_pairs)
                cmd = [sys.executable,'-u',str(engine),'--root',str(root),'--sample',lib['sample'],
                    '--host-fasta',str(host),'--threads',str(threads),'--chunk-pairs',str(chunk),
                    '--r1',*map(str,lib['r1']),'--r2',*map(str,lib['r2'])]
                if resume: cmd.append('--resume')
                if args.allow_index_build: cmd.append('--allow-index-build')
                for name in ('condition','biological_replicate'):
                    if cfg.get(name): cmd.extend(['--'+name.replace('_','-'),str(cfg[name])])
                expected_config = {**THRESHOLDS,'sample':lib['sample'],'side':lib['side'],
                    'analysis_version':version,'threads':threads,'chunk_pairs':chunk,
                    '_source_sha256':{name:sha256_file(engine.parent/name) for name in SOURCE_FILES}}
                print(f"\n[{i}/{len(libraries)}] {lib['sample']} — {s['mode']}",flush=True)
                checkpoint(outdir,manifest,statuses,'running',started)
                rc = run_child(cmd,Path(s['log_file']))
                if rc!=0: raise ValueError(f'Analysis exited with code {rc}; see {s["log_file"]}')
                s['validated'] = validate_result(lib['sample'],lib['side'],result_dir,expected_config,lib)
                s['status'] = 'complete'
            except KeyboardInterrupt:
                s.update(status='interrupted',error='Interrupted by user'); interrupted=True
            except (ValueError,OSError,KeyError,TypeError) as e:
                s.update(status='failed',error=str(e))
                print('FAILED '+lib['sample']+': '+str(e),file=sys.stderr,flush=True)
            finally:
                s['elapsed_seconds'] = round(time.perf_counter()-begin,3)
                checkpoint(outdir,manifest,statuses,'interrupted' if interrupted else 'running',started)
            if interrupted: break
        try:
            summary = aggregate_results(libraries,statuses,outdir)
        except (ValueError,OSError,KeyError,TypeError) as e:
            checkpoint(outdir,manifest,statuses,'aggregation_failed',started,{'error':str(e)})
            raise
        state = 'interrupted' if interrupted else summary['status']
        checkpoint(outdir,manifest,statuses,state,started,summary if not interrupted else {**summary,'status':'interrupted'})
        print(f"\nBATCH {state.upper()}: {summary['completed_libraries']}/{len(libraries)} libraries complete. {outdir}",flush=True)
        if state!='complete': print('Only partial tables are available. Check library_status.csv, then rerun the same command.',flush=True)
        return {**summary,'status':state}


def main(argv=None):
    try:
        result = run_batch(build_parser().parse_args(argv))
        return 0 if result['status'] in ('complete','plan_only') else 130 if result['status']=='interrupted' else 1
    except KeyboardInterrupt:
        print('Interrupted; rerun the same batch command.',file=sys.stderr); return 130
    except (ValueError,OSError,KeyError,TypeError) as e:
        print('ERROR: '+str(e),file=sys.stderr); return 1


if __name__=='__main__':
    raise SystemExit(main())
