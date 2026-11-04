#!/usr/bin/env python3
"""Read-only, auditable preparation of library-level RAISING figure data.

Sequence groups are candidate patterns, not independent integrations. Original
annotations are retained; insert-containing junction MH is unavailable. This
exporter does not infer experimental metadata from sample names. The default
bin limits describe the bundled 800-bp study profile, not biological maxima.
"""
import argparse
import collections
import csv
import gzip
import hashlib
import json
from pathlib import Path
import statistics

VERSION = '1.0.0'
STATUS_KEYS = ('template_supported','template_ambiguous','alignment_ambiguous',
               'junction_variant_review','no_qualifying_template','no_nearby_template')
GROUP_NUMERIC = ('supported_pairs','review_pairs','without_prior_dd_or_complex_host_pairs',
                 'prior_dd_without_complex_host_pairs','complex_host_pairs')
ANNOTATION_NUMERIC = ('host_boundary0','host_boundary_min0','host_boundary_max0',
    'donor_deletion_min_bp','donor_deletion_max_bp','mh_length_bp','insertion_length_bp',
    'pair_support','simple_pair_support','prior_dd_pair_support','complex_host_pair_support',
    'is_group_representative','supported_pair_count_for_identical_sequence_group',
    'supported_pair_count_once_per_group','review_pair_count_for_identical_sequence_group',
    'template_reference_mh_length_bp','template_reference_donor_deletion_min_bp',
    'template_reference_donor_deletion_max_bp')
ANNOTATION_KEEP = ('sample','library_side','condition','biological_replicate','event_id',
    'chrom','host_strand','equivalent_template_group','reference_compartment','mh_sequence',
    'insertion_sequence','mh_for_analysis_bp','mh_assessment',*ANNOTATION_NUMERIC)
COMPARISON_KEEP = ('sample','template_id','reference_template_id','relationship',
    'has_junction_proximal_difference','differences_only_donor_internal_substitutions',
    'alignment_ambiguous','alignment_edit_bases','substitutions','indel_bases','gap_events',
    'operations_json')
SUBSETS = ('nuclear_all','nuclear_with_simple','nuclear_only_simple')
TABLE_NAMES = ('all_identical_sequence_groups.tsv', 'all_junction_reassessment.tsv',
               'all_template_comparisons.tsv')

def require(condition, message):
    """Use explicit checks so validation cannot disappear under python -O."""
    if not condition:
        raise ValueError(message)

def input_path(directory, name):
    candidates = [directory / name, directory / (name + '.gz')]
    present = [p for p in candidates if p.is_file()]
    require(len(present) == 1,
            f'Expected exactly one of {name} or {name}.gz; found {len(present)}')
    return present[0]

def rows(path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', encoding='utf-8-sig', newline='') as f:
        yield from csv.DictReader(f,delimiter='\t')

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()

def metadata(sample, side, values):
    return {'sample': sample, 'library_side': side, 'series_code': '',
            'replicate_label_from_name': '',
            'condition': values.get('condition', ''),
            'biological_replicate': values.get('biological_replicate', ''),
            'metadata_status': ('metadata_provided' if any(values.values())
                                else 'metadata_not_provided')}

def merge_metadata(source, libraries, values):
    sample = source['sample']
    require(sample in libraries, f'Unknown sample in table or metadata: {sample}')
    side = source.get('library_side', '')
    require(not side or side == libraries[sample]['library_side'],
            f'Library side conflict for {sample}')
    for key in ('condition', 'biological_replicate'):
        value = source.get(key, '')
        value = '' if value is None else str(value)
        previous = values[sample].get(key, '')
        require(not (previous and value and previous != value),
                f'Conflicting {key} for {sample}')
        if value:
            values[sample][key] = value

def fraction(n,d):return n/d if d else None

def write_csv(path,items,fields=None):
    fields = list(items[0]) if items else list(fields or ())
    require(bool(fields), f'No schema for empty output {path.name}')
    with path.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fields);w.writeheader();w.writerows(items)

def selected_basis(row,subset):
    if row['reference_compartment']!='nuclear_or_unplaced_reference':return 0
    if subset=='nuclear_all':return row['supported_pairs']
    if subset=='nuclear_only_simple' and (row['prior_dd_without_complex_host_pairs'] or row['complex_host_pairs']):return 0
    return row['without_prior_dd_or_complex_host_pairs']

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-dir',type=Path,required=True,
                   help='Complete junction batch directory; TSV or TSV.GZ tables accepted')
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--metadata',type=Path,
                   help='Optional CSV with sample, library_side, condition, biological_replicate; conflicting values fail')
    p.add_argument('--max-donor-loss-3p',type=int,default=39,
                   help='Expected export bound for the analysis profile (default:39); does not rerun junction calling')
    p.add_argument('--max-donor-loss-5p',type=int,default=45,
                   help='Expected export bound for the analysis profile (default:45); does not rerun junction calling')
    p.add_argument('--max-insertion',type=int,default=20,
                   help='Expected insertion search/export bound (default:20); does not rerun junction calling')
    a=p.parse_args(argv)
    require(min(a.max_donor_loss_3p, a.max_donor_loss_5p, a.max_insertion) >= 0,
            'Expected analysis bounds must be nonnegative')
    paths={n:input_path(a.input_dir,n) for n in TABLE_NAMES}
    paths['batch_summary.json']=a.input_dir/'batch_summary.json'
    if a.metadata:
        paths['metadata.csv']=a.metadata
    inputs={n:{'file_name':path.name,'bytes':path.stat().st_size,'sha256':digest(path)}
            for n,path in paths.items()}
    summary=json.loads((a.input_dir/'batch_summary.json').read_text())
    require(summary['status']=='complete', 'Batch must have status=complete')
    libraries={x['sample']:x for x in summary['libraries']}
    require(len(libraries) == len(summary['libraries']) > 0, 'Duplicate or empty library inventory')
    require(len(libraries)==summary['completed_library_count']==summary['selected_library_count'],
            'Completed/selected library counts do not reconcile')
    require(all(x['status']=='complete' for x in libraries.values()), 'Incomplete library in batch')
    require(all(x['library_side'] in ('3p','5p') for x in libraries.values()), 'Invalid library side')
    metadata_values = {sample:{} for sample in libraries}
    for lib in libraries.values():
        merge_metadata(lib, libraries, metadata_values)
    if a.metadata:
        with a.metadata.open(encoding='utf-8-sig', newline='') as f:
            reader=csv.DictReader(f)
            require({'sample','library_side','condition','biological_replicate'} <= set(reader.fieldnames or ()),
                    'Metadata CSV needs sample, library_side, condition, biological_replicate columns')
            seen=set()
            for row in reader:
                require(row['sample'] not in seen, 'Duplicate metadata sample: '+row['sample'])
                seen.add(row['sample'])
                merge_metadata(row, libraries, metadata_values)
    groups={}
    for r in rows(paths['all_identical_sequence_groups.tsv']):
        merge_metadata(r, libraries, metadata_values)
        key=(r['sample'],r['equivalent_template_group']);require(key not in groups, 'Validation failed: key not in groups')
        for k in GROUP_NUMERIC:r[k]=int(r[k])
        require(all(r[k] >= 0 for k in GROUP_NUMERIC), 'Negative sequence-group support count')
        require(r['supported_pairs']==sum(r[k] for k in GROUP_NUMERIC[2:]), "Validation failed: r['supported_pairs']==sum(r[k] for k in GROUP_NUMERIC[2:])")
        require(r['template_ids']==r['representative_template_id'], 'This export requires single-member groups; retain competing annotations if a future run differs.')
        groups[key]=r
    annotations={};by_sample=collections.defaultdict(list)
    for source in rows(paths['all_junction_reassessment.tsv']):
        merge_metadata(source, libraries, metadata_values)
        r={k:source[k] for k in ANNOTATION_KEEP}
        for k in ANNOTATION_NUMERIC:r[k]=int(r[k])
        require(all(r[k] >= 0 for k in ANNOTATION_NUMERIC), 'Negative coordinate, length, or count')
        key=(r['sample'],r['equivalent_template_group']);g=groups[key]
        require(r['event_id']==g['representative_template_id'] and r['is_group_representative']==1, "Validation failed: r['event_id']==g['representative_template_id'] and r['is_group_representative']==1")
        require(r['reference_compartment']==g['reference_compartment'], "Validation failed: r['reference_compartment']==g['reference_compartment']")
        require(r['reference_compartment'] in ('nuclear_or_unplaced_reference','mitochondrial_reference'),
                'Unsupported reference compartment')
        require(r['supported_pair_count_once_per_group']==g['supported_pairs'], "Validation failed: r['supported_pair_count_once_per_group']==g['supported_pairs']")
        require(r['supported_pair_count_for_identical_sequence_group']==g['supported_pairs'], "Validation failed: r['supported_pair_count_for_identical_sequence_group']==g['supported_pairs']")
        require(r['review_pair_count_for_identical_sequence_group']==g['review_pairs'], "Validation failed: r['review_pair_count_for_identical_sequence_group']==g['review_pairs']")
        require(r['donor_deletion_max_bp']-r['donor_deletion_min_bp']==r['mh_length_bp'], "Validation failed: r['donor_deletion_max_bp']-r['donor_deletion_min_bp']==r['mh_length_bp']")
        require(len(r['mh_sequence'])==r['mh_length_bp'], "Validation failed: len(r['mh_sequence'])==r['mh_length_bp']")
        require(len(r['insertion_sequence'])==r['insertion_length_bp'], "Validation failed: len(r['insertion_sequence'])==r['insertion_length_bp']")
        require(r['template_reference_mh_length_bp']==r['mh_length_bp'], "Validation failed: r['template_reference_mh_length_bp']==r['mh_length_bp']")
        require(r['template_reference_donor_deletion_min_bp']==r['donor_deletion_min_bp'], "Validation failed: r['template_reference_donor_deletion_min_bp']==r['donor_deletion_min_bp']")
        require(r['template_reference_donor_deletion_max_bp']==r['donor_deletion_max_bp'], "Validation failed: r['template_reference_donor_deletion_max_bp']==r['donor_deletion_max_bp']")
        require(r['host_boundary_max0']-r['host_boundary_min0']==r['mh_length_bp'], "Validation failed: r['host_boundary_max0']-r['host_boundary_min0']==r['mh_length_bp']")
        if r['insertion_length_bp']:
            require(r['mh_for_analysis_bp']=='' and r['mh_assessment']=='not_assessed_insertion_interfaces', "Validation failed: r['mh_for_analysis_bp']=='' and r['mh_assessment']=='not_assessed_insertion_interfaces'")
            r['mh_for_analysis_bp']=None
        else:
            require(int(r['mh_for_analysis_bp'])==r['mh_length_bp'], "Validation failed: int(r['mh_for_analysis_bp'])==r['mh_length_bp']")
            r['mh_for_analysis_bp']=int(r['mh_for_analysis_bp'])
        r.update({k:g[k] for k in GROUP_NUMERIC})
        ekey=(r['sample'],r['event_id']);require(ekey not in annotations, 'Validation failed: ekey not in annotations')
        annotations[ekey]=r;by_sample[r['sample']].append(r)
    comparisons={};relationship_counts=collections.Counter()
    for source in rows(paths['all_template_comparisons.tsv']):
        merge_metadata(source, libraries, metadata_values)
        r={k:source[k] for k in COMPARISON_KEEP}
        key=(r['sample'],r['template_id']);require(key not in comparisons and key in annotations, 'Validation failed: key not in comparisons and key in annotations')
        for k in ('alignment_edit_bases','substitutions','indel_bases','gap_events'):r[k]=int(r[k])
        for k in ('has_junction_proximal_difference','differences_only_donor_internal_substitutions','alignment_ambiguous'):
            require(r[k] in ('True','False'), "Validation failed: r[k] in ('True','False')");r[k]=(r[k]=='True')
        comparisons[key]=r;relationship_counts[r['relationship']]+=1
        annotations[key]['comparison_relationship']=r['relationship']
        annotations[key]['alternative_representation_flag']=r['relationship'] in ('substitution_only_alternative','small_indel_alternative')
    require(set(comparisons)==set(annotations), 'Validation failed: set(comparisons)==set(annotations)')
    require(len(groups)==len(annotations), 'Group/annotation counts differ')
    for name,count in zip(TABLE_NAMES,(len(groups),len(annotations),len(comparisons))):
        recorded=summary['merged_outputs'].get(name+'.gz',summary['merged_outputs'].get(name))
        require(recorded is not None and recorded['rows']==count, 'Merged row count mismatch: '+name)
        if paths[name].suffix=='.gz' and recorded.get('sha256'):
            require(inputs[name]['sha256']==recorded['sha256'], 'Compressed source hash mismatch: '+name)
    require(set(by_sample) <= set(libraries), 'Annotations contain an unlisted sample')
    for r in annotations.values():
        r.update(metadata_values[r['sample']])
    qc=[];metrics=[];mh_hist=[];del_hist=[];ins_hist=[];alternative_review=[]
    max_mh=max((r['mh_length_bp'] for r in annotations.values()),default=0)
    loss_bounds={'3p':a.max_donor_loss_3p,'5p':a.max_donor_loss_5p}
    for r in annotations.values():
        require(r['donor_deletion_max_bp'] <= loss_bounds[r['library_side']],
                'Donor loss exceeds the specified analysis-profile bound')
        require(r['insertion_length_bp'] <= a.max_insertion,
                'Insertion length exceeds the specified analysis-profile bound')
    for sample,lib in sorted(libraries.items()):
        rs=by_sample[sample];side=lib['library_side'];meta=metadata(sample,side,metadata_values[sample])
        require(all(isinstance(lib[k], int) and lib[k] >= 0 for k in
                    (*STATUS_KEYS, 'total_input_pairs', 'host_supported_pairs',
                     'original_candidate_patterns', 'identical_sequence_groups',
                     'supported_identical_sequence_groups')),
                'Batch counts must be nonnegative integers: '+sample)
        require(lib['template_supported'] <= lib['host_supported_pairs'] <= lib['total_input_pairs'],
                'Read-pair denominators do not reconcile: '+sample)
        require(len(rs)==lib['original_candidate_patterns']==lib['identical_sequence_groups'], "Validation failed: len(rs)==lib['original_candidate_patterns']==lib['identical_sequence_groups']")
        supported=[r for r in rs if r['supported_pairs']>0]
        nuclear=[r for r in supported if r['reference_compartment']=='nuclear_or_unplaced_reference']
        mt=[r for r in supported if r['reference_compartment']=='mitochondrial_reference']
        require(sum(r['supported_pairs'] for r in rs)==lib['template_supported'], "Validation failed: sum(r['supported_pairs'] for r in rs)==lib['template_supported']")
        require(len(supported)==lib['supported_identical_sequence_groups'], "Validation failed: len(supported)==lib['supported_identical_sequence_groups']")
        require(sum(lib[k] for k in STATUS_KEYS)==lib['host_supported_pairs'], "Validation failed: sum(lib[k] for k in STATUS_KEYS)==lib['host_supported_pairs']")
        for comp,rr in [('nuclear_or_unplaced_reference',nuclear),('mitochondrial_reference',mt)]:
            require(sum(r['supported_pairs'] for r in rr)==lib['supported_pairs_by_reference_compartment'].get(comp,0), "Validation failed: sum(r['supported_pairs'] for r in rr)==lib['supported_pairs_by_reference_compartment'].get(comp,0)")
        qc.append({**meta,'raw_pairs':lib['total_input_pairs'],'host_supported_pairs':lib['host_supported_pairs'],
                   **{k:lib[k] for k in STATUS_KEYS},'candidate_patterns_total':len(rs),
                   'supported_patterns_total':len(supported),'nuclear_supported_patterns':len(nuclear),
                   'MT_supported_patterns':len(mt),'nuclear_supported_pairs':sum(r['supported_pairs'] for r in nuclear),
                   'MT_supported_pairs':sum(r['supported_pairs'] for r in mt),
                   'singleton_patterns':sum(r['supported_pairs']==1 for r in supported),
                   'patterns_with_insertion_candidate':sum(r['insertion_length_bp']>0 for r in supported),
                   'patterns_with_review_pairs':sum(r['review_pairs']>0 for r in supported),
                   'alternative_representation_patterns':sum(r['alternative_representation_flag'] for r in supported),
                   'host_mapping_fraction_raw':fraction(lib['host_supported_pairs'],lib['total_input_pairs']),
                   'junction_supported_fraction_host':fraction(lib['template_supported'],lib['host_supported_pairs']),
                   'singleton_fraction_patterns':fraction(sum(r['supported_pairs']==1 for r in supported),len(supported))})
        for subset in SUBSETS:
            for threshold in (1,2):
                selected=[r for r in rs if selected_basis(r,subset)>=threshold]
                direct=[r for r in selected if r['insertion_length_bp']==0]
                ins=[r for r in selected if r['insertion_length_bp']>0]
                md={**meta,'subset':subset,'min_support':threshold,
                    'support_basis':'all_supported_pairs' if subset=='nuclear_all' else 'without_prior_dd_or_complex_host_pairs'}
                record={**md,'N_patterns':len(selected),'N_direct':len(direct),'N_insertion_candidate':len(ins),
                    'N_direct_mh0':sum(r['mh_length_bp']==0 for r in direct),
                    'N_direct_mh1':sum(r['mh_length_bp']==1 for r in direct),
                    'N_direct_mh2plus':sum(r['mh_length_bp']>=2 for r in direct),
                    'N_direct_deletion_min_gt0':sum(r['donor_deletion_min_bp']>0 for r in direct),
                    'N_direct_deletion_zero_exact':sum(r['donor_deletion_min_bp']==r['donor_deletion_max_bp']==0 for r in direct),
                    'N_direct_deletion_ambiguous_zero':sum(r['donor_deletion_min_bp']==0 and r['donor_deletion_max_bp']>0 for r in direct),
                    'N_reviewflagged_patterns':sum(r['review_pairs']>0 for r in selected),
                    'N_alternative_representation':sum(r['alternative_representation_flag'] for r in selected),
                    'selectedsupport_pairs':sum(selected_basis(r,subset) for r in selected),
                    'totalgroup_support_pairs':sum(r['supported_pairs'] for r in selected),
                    'singleton_totalpairs_count':sum(r['supported_pairs']==1 for r in selected),
                    'direct_deletion_min_median_bp':statistics.median(r['donor_deletion_min_bp'] for r in direct) if direct else None,
                    'direct_deletion_max_median_bp':statistics.median(r['donor_deletion_max_bp'] for r in direct) if direct else None,
                    'direct_deletion_min_sum_bp':sum(r['donor_deletion_min_bp'] for r in direct),
                    'direct_deletion_max_sum_bp':sum(r['donor_deletion_max_bp'] for r in direct)}
                require(record['N_patterns']==record['N_direct']+record['N_insertion_candidate'], "Validation failed: record['N_patterns']==record['N_direct']+record['N_insertion_candidate']")
                require(record['N_direct']==sum(record[k] for k in ('N_direct_mh0','N_direct_mh1','N_direct_mh2plus')), "Validation failed: record['N_direct']==sum(record[k] for k in ('N_direct_mh0','N_direct_mh1','N_direct_mh2plus'))")
                require(record['N_direct']==sum(record[k] for k in ('N_direct_deletion_min_gt0','N_direct_deletion_zero_exact','N_direct_deletion_ambiguous_zero')), "Validation failed: record['N_direct']==sum(record[k] for k in ('N_direct_deletion_min_gt0','N_direct_deletion_zero_exact','N_direct_deletion_ambiguous_zero'))")
                metrics.append(record)
                for bp in range(max_mh+1):
                    binrows=[r for r in direct if r['mh_length_bp']==bp]
                    mh_hist.append({**md,'mh_bp':bp,'pattern_count':len(binrows),'direct_pattern_denominator':len(direct),
                        'selectedsupport_pairs_in_bin':sum(selected_basis(r,subset) for r in binrows),
                        'direct_selectedsupport_pair_denominator':sum(selected_basis(r,subset) for r in direct),
                        'fraction_patterns':fraction(len(binrows),len(direct))})
                max_bp=loss_bounds[side]
                require(all(r['donor_deletion_max_bp']<=max_bp for r in direct), "Validation failed: all(r['donor_deletion_max_bp']<=max_bp for r in direct)")
                for bp in range(max_bp+1):
                    mn=[r for r in direct if r['donor_deletion_min_bp']==bp]
                    mx=[r for r in direct if r['donor_deletion_max_bp']==bp]
                    del_hist.append({**md,'deletion_bp':bp,'count_min_endpoint':len(mn),'count_max_endpoint':len(mx),
                        'direct_pattern_denominator':len(direct),'selectedsupport_pairs_min_endpoint':sum(selected_basis(r,subset) for r in mn),
                        'selectedsupport_pairs_max_endpoint':sum(selected_basis(r,subset) for r in mx),
                        'direct_selectedsupport_pair_denominator':sum(selected_basis(r,subset) for r in direct),
                        'fraction_min_endpoint':fraction(len(mn),len(direct)),
                        'fraction_max_endpoint':fraction(len(mx),len(direct))})
                for bp in range(1,a.max_insertion+1):
                    br=[r for r in ins if r['insertion_length_bp']==bp]
                    ins_hist.append({**md,'insertion_bp':bp,'pattern_count':len(br),'insertion_pattern_denominator':len(ins),
                        'selectedsupport_pairs_in_bin':sum(selected_basis(r,subset) for r in br),
                        'insertion_selectedsupport_pair_denominator':sum(selected_basis(r,subset) for r in ins),
                        'fraction_patterns':fraction(len(br),len(ins))})
        for r in rs:
            comp=comparisons[(sample,r['event_id'])]
            if comp['relationship']=='strongest_template':continue
            alternative_review.append({**meta,'event_id':r['event_id'],
                'equivalent_template_group':r['equivalent_template_group'],
                'reference_compartment':r['reference_compartment'],'chrom':r['chrom'],'host_strand':r['host_strand'],
                'host_boundary0':r['host_boundary0'],'donor_deletion_min_bp':r['donor_deletion_min_bp'],
                'donor_deletion_max_bp':r['donor_deletion_max_bp'],'insertion_length_bp':r['insertion_length_bp'],
                'insertion_sequence':r['insertion_sequence'],'mh_for_analysis_bp':r['mh_for_analysis_bp'],
                **{k:r[k] for k in GROUP_NUMERIC},**{k:v for k,v in comp.items() if k not in ('sample','template_id')},
                'original_annotation_preserved':True})
    figure_values = []
    for record in metrics:
        figure_values.append({**record,
            'MH0_fraction':fraction(record['N_direct_mh0'],record['N_direct']),
            'MH1_fraction':fraction(record['N_direct_mh1'],record['N_direct']),
            'MH_ge2_fraction':fraction(record['N_direct_mh2plus'],record['N_direct']),
            'donor_loss_min_gt0_fraction':fraction(record['N_direct_deletion_min_gt0'],record['N_direct']),
            'insertion_candidate_fraction':fraction(record['N_insertion_candidate'],record['N_patterns'])})
    outputs={'figure_values':figure_values,'metrics':metrics,'qc':qc,'mh_distribution':mh_hist,'deletion_distribution':del_hist,
             'insertion_distribution':ins_hist,'alternative_review':alternative_review}
    for name,source in paths.items():
        require(digest(source)==inputs[name]['sha256'], 'Source changed during processing: '+name)
    a.output_dir.mkdir(parents=True,exist_ok=True)
    empty_schemas={
        'alternative_review': [*metadata('', '', {}).keys(), 'event_id', 'equivalent_template_group',
            'reference_compartment','chrom','host_strand','host_boundary0','donor_deletion_min_bp',
            'donor_deletion_max_bp','insertion_length_bp','insertion_sequence','mh_for_analysis_bp',
            *GROUP_NUMERIC, *[k for k in COMPARISON_KEEP if k not in ('sample','template_id')],
            'original_annotation_preserved'],
        'insertion_distribution': [*metadata('', '', {}).keys(),'subset','min_support','support_basis',
            'insertion_bp','pattern_count','insertion_pattern_denominator','selectedsupport_pairs_in_bin',
            'insertion_selectedsupport_pair_denominator','fraction_patterns']}
    for name,items in outputs.items():
        write_csv(a.output_dir/(name+'.csv'),items,empty_schemas.get(name))
    # A lightweight candidate provenance export enables independent filtering.
    write_csv(a.output_dir/'candidate_annotations.csv',list(annotations.values()),
              [*ANNOTATION_KEEP,*GROUP_NUMERIC,'comparison_relationship','alternative_representation_flag'])
    payload={'analysis_version':VERSION,'sources':inputs,'batch_status':summary['status'],
        'counting_unit':'One supported candidate sequence group per sample; not an independently validated integration event.',
        'metadata_policy':'Condition and biological replicate are retained from consistent explicit input metadata. Optional metadata CSV may fill blanks; conflicts fail. Sample names are not parsed. Compatibility columns series_code and replicate_label_from_name remain blank.',
        'subset_definitions':{
            'nuclear_all':'Nuclear/unplaced reference; threshold applies to total supported pairs.',
            'nuclear_with_simple':'Nuclear/unplaced reference; threshold applies to simple-support pairs. Groups may also contain prior-DD or complex support.',
            'nuclear_only_simple':'Nuclear/unplaced reference; no supported prior-DD or complex pairs; threshold applies to simple-support pairs.'},
        'annotation_policy':'Use inherited template descriptions. Insertion interfaces have MH NA. Substitution/small-indel alternatives remain flagged, not corrected or removed.',
        'limits':{'donor_deletion_nominal_max_3p':a.max_donor_loss_3p,
                  'donor_deletion_nominal_max_5p':a.max_donor_loss_5p,'insertion_max':a.max_insertion,
                  'limits_are_expected_analysis_bounds_not_biological_maxima':True,
                  'MH_observed_max':max_mh,'MH_observed_max_is_not_detection_guarantee':True},
        'checks':{'all_libraries_complete':True,'all_joins_match':True,'all_groups_single_member':True,
                  'per_library_denominators_reconcile':True,'sequence_and_boundary_invariants':True,
                  'all_metric_categories_reconcile':True,'explicit_metadata_consistent':True,
                  'original_source_files_modified':False,
                  'compressed_hash_verification':'Checked where compressed source bytes and expected SHA-256 are available. Uncompressed inputs have their own SHA-256 recorded; they cannot be compared to a compressed-file hash.'},
        'totals':{'libraries':len(libraries),'candidate_groups':len(groups),'supported_pairs':sum(q['template_supported'] for q in qc),
                  'supported_patterns':sum(q['supported_patterns_total'] for q in qc),'relationship_counts':dict(relationship_counts)},
        'output_rows':{name:len(items) for name,items in outputs.items()},
        'candidate_annotation_rows':len(annotations)}
    output=a.output_dir/'analysis_summary.json'
    output.write_text(json.dumps(payload,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'json':str(output),'csv_dir':str(a.output_dir),'rows':payload['output_rows'],
                      'totals':payload['totals'],'observed_max_mh':max_mh},indent=2))
    return payload

if __name__=='__main__':main()
