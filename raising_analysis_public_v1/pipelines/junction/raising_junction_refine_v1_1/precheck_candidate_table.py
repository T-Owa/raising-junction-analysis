#!/usr/bin/env python3
"""Review an existing junction table without claiming raw-read realignment."""
import argparse
import collections
import csv
import hashlib
import json
from pathlib import Path
import sys

import refine_templates as rt


def write_table(path, rows, fallback):
    fields = list(dict.fromkeys(key for row in rows for key in row)) or fallback
    with path.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fields, delimiter='\t')
        w.writeheader()
        for r in rows:
            w.writerow({k: json.dumps(v, sort_keys=True, ensure_ascii=False) if isinstance(v, (list, dict)) else v
                        for k, v in r.items()})


def run(junctions, summary, donor_path, output):
    junctions, summary, donor_path, output = map(Path, (junctions, summary, donor_path, output))
    manifest = json.loads(summary.read_text(encoding='utf-8'))
    digest = hashlib.sha256(junctions.read_bytes()).hexdigest()
    if manifest.get('status') != 'complete' or digest != manifest['output_sha256']['junctions.tsv']:
        raise ValueError('Pilot must be complete and junction table must match its SHA-256')
    donor = ''.join(line.strip() for line in donor_path.read_text().splitlines() if not line.startswith('>')).upper()
    if hashlib.sha256(donor.encode()).hexdigest() != manifest['reference']['donor_sha256']:
        raise ValueError('Donor sequence does not match pilot reference')
    with junctions.open(newline='') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    clusters = rt.cluster_candidates(rows, max_span=10)
    alternatives = rt.build_donor_only_review(rows, donor, manifest['library_side'])
    group_rows = []
    for cluster in clusters:
        members = cluster['rows']
        group_rows.append({k: v for k, v in cluster.items() if k != 'rows'} | {
            'candidate_patterns': len(members),
            'pilot_supporting_pairs': sum(int(r['pair_support']) for r in members),
            'event_ids': [r['event_id'] for r in members],
            'interpretation': 'Review group only; not a validated insertion event'})
    simple = [r for r in rows if int(r['simple_pair_support']) > 0]
    report = {
        'analysis_type': 'junction_candidate_table_precheck_not_raw_realignment',
        'sample': manifest['sample'], 'library_side': manifest['library_side'],
        'input_junction_sha256': digest,
        'pilot_stats': manifest['stats'],
        'max_group_span_bp': 10,
        'all_candidate_patterns': len(rows),
        'review_groups_all': len(clusters),
        'multipattern_review_groups_all': sum(len(c['rows']) > 1 for c in clusters),
        'simple_supported_patterns': len(simple),
        'simple_single_pair_patterns': sum(int(r['simple_pair_support']) == 1 for r in simple),
        'simple_at_least_two_pair_patterns': sum(int(r['simple_pair_support']) >= 2 for r in simple),
        'simple_exact_host_boundary_keys': len({(r['chrom'], r['host_strand'], r['host_boundary0']) for r in simple}),
        'mt_patterns_all': sum(r['chrom'] in ('MT','chrM','M','chrMT') for r in rows),
        'mt_patterns_simple': sum(r['chrom'] in ('MT','chrM','M','chrMT') for r in simple),
        'donor_only_comparison_rows': len(alternatives),
        'comparison_labels': dict(collections.Counter(r.get('relationship', r.get('classification', r.get('comparison_type', 'see_table')))
                                                       for r in alternatives)),
        'limitations': [
            'No full raw FASTQ, qualities, or human reference realignment was performed in this table-only precheck.',
            'Donor-only comparisons are restricted to exactly the same canonical human boundary; near-boundary comparisons need host sequence.',
            'Small sequence differences do not establish PCR/sequencing error and do not justify merging genuine terminal deletions.',
            'Review groups and candidate patterns are not independent molecule, cell, or validated integration counts.',
        ],
    }
    output.mkdir(parents=True, exist_ok=True)
    write_table(output/'boundary_review_groups.tsv', group_rows, ['cluster_id'])
    write_table(output/'same_boundary_donor_alternatives.tsv', alternatives, ['event_id'])
    (output/'table_precheck_summary.json').write_text(json.dumps(report, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps({k: v for k,v in report.items() if k not in ('limitations','pilot_stats')}, indent=2))
    return report


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--junctions', required=True)
    p.add_argument('--summary', required=True)
    p.add_argument('--donor', required=True)
    p.add_argument('--output-dir', required=True)
    a = p.parse_args()
    try:
        run(a.junctions, a.summary, a.donor, a.output_dir)
    except (ValueError, KeyError, OSError) as e:
        print('ERROR: '+str(e), file=sys.stderr)
        raise SystemExit(2)
