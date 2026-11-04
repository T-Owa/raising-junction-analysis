"""Auditable four-class exports for the frozen RAISING screen (v1.0.0)."""
import collections
import csv
import gzip
import hashlib
import json
import re
from pathlib import Path

import raising_threeway_classify as core
import raising_concatemer_screen as ds

VERSION = '1.0.0'
CLASSES = (
    'donor_rearrangement_supported', 'host_mappable_other',
    'short_donor_primer_candidate', 'other_unresolved',
)
LABELS = {
    CLASSES[0]: 'Donor junction/rearrangement supported (concatemer candidate)',
    CLASSES[1]: 'Other host-mapping supported',
    CLASSES[2]: 'Short donor/primer product candidate',
    CLASSES[3]: 'Other/unresolved',
}
FILES = ('pair_classification.tsv', 'read_evidence.tsv',
         'donor_junction_evidence.tsv', 'host_alignment_evidence.tsv')


def short_product_pattern(raw_sequence, donor, side='3p'):
    """R1 pattern only; not a monomer, donor physical end or deletion call.

    Exact >=15-nt prefix at P3 start (allow 0-2 nt missing at the read start),
    followed by the first maximal A run of >=15. No mismatches in prefix.
    A-only patterns are retained as candidates, with confirmation separate.
    """
    seq = raw_sequence.upper()
    oriented = donor if side == '3p' else ds.revcomp(donor)
    primer = ds.P3[side]
    if oriented.count(primer) != 1:
        raise ValueError('P3 is not unique in the oriented donor')
    start = oriented.index(primer)
    result = {'short_pattern': False, 'short_confidence': '',
              'prefix_bp': 0, 'prefix_donor_oriented_start_1based': '',
              'bases_beyond_P3': 0, 'polyA_bp': 0,
              'adapter_confirmed': False, 'adapter_exact': False,
              'primer_only_prefix': False}
    run = re.search('A{15,}', seq)
    if run is None or run.start() < 15:
        return result
    prefix = seq[:run.start()]
    starts = [s for s in range(start, start+3)
              if oriented[s:s+len(prefix)] == prefix]
    if not starts:
        return result
    s = starts[0]
    exact = seq[run.end():].startswith(core.R1_ADAPTER)
    confirmed = any(
        end-run.start() >= 15 and
        core.hamming(seq[end:end+len(core.R1_ADAPTER)], core.R1_ADAPTER) <= 2
        for end in (run.end(), run.end()-1)
    )
    beyond = max(0, s+len(prefix)-(start+len(primer)))
    result.update(short_pattern=True,
                  short_confidence='adapter_confirmed' if confirmed else 'polyA_only_candidate',
                  prefix_bp=len(prefix), prefix_donor_oriented_start_1based=s+1,
                  bases_beyond_P3=beyond, polyA_bp=run.end()-run.start(),
                  adapter_confirmed=confirmed, adapter_exact=exact,
                  primer_only_prefix=beyond==0)
    return result


def donor_join_subtype(rec):
    """Selected donor anchors only; host interval containment is a candidate flag."""
    if not rec['dd']:
        return '', 0
    d = rec['donor_details']
    gap = int(d['unaligned_gap_bp'])
    if gap <= 0:
        return 'abutting_or_overlapping_donor_anchors', 0
    lo = int(d['left_read_end_1based'])+1
    hi = int(d['right_read_start_1based'])-1
    best_overlap = 0
    contained = False
    for h in rec['host_hits']:
        if not h.get('accepted'):
            continue
        hs, he = h.get('raw_query_start_1based'), h.get('raw_query_end_1based')
        if not hs or not he:
            continue
        hs, he = int(hs), int(he)
        best_overlap = max(best_overlap, max(0, min(hi, he)-max(lo, hs)+1))
        contained |= lo <= hs <= he <= hi
    return ('host_interval_within_donor_gap_candidate' if contained
            else 'intervening_sequence_between_donor_anchors'), best_overlap


def classify_pair(ends, pattern):
    old, dd, host, reasons = core.assign_pair(ends)
    category = (CLASSES[0] if dd else CLASSES[1] if host else
                CLASSES[2] if pattern['short_pattern'] else CLASSES[3])
    subtypes = {donor_join_subtype(r)[0] for r in ends if r['dd']}
    return {'classification': category, 'original_threeway_class': old,
            'donor_rearrangement_supported': int(dd),
            'host_mapping_supported': int(host), 'both_supported': int(dd and host),
            'same_read_both_supported': int(any(r['dd'] and r['host'] for r in ends)),
            'abutting_or_overlapping_donor_anchors': int('abutting_or_overlapping_donor_anchors' in subtypes),
            'donor_gap_with_host_interval_candidate': int('host_interval_within_donor_gap_candidate' in subtypes),
            'donor_gap_other_candidate': int('intervening_sequence_between_donor_anchors' in subtypes),
            'junction_types': ';'.join(sorted({r['donor_details'].get('junction_type','') for r in ends if r['dd']})),
            'R1_short_pattern': int(pattern['short_pattern']),
            'R1_short_confidence': pattern['short_confidence'],
            'R1_prefix_bp': pattern['prefix_bp'],
            'R1_prefix_donor_oriented_start_1based': pattern['prefix_donor_oriented_start_1based'],
            'R1_bases_beyond_P3': pattern['bases_beyond_P3'],
            'R1_polyA_bp': pattern['polyA_bp'],
            'R1_adapter_confirmed': int(pattern['adapter_confirmed']),
            'R1_adapter_exact': int(pattern['adapter_exact']),
            'R1_primer_only_prefix': int(pattern['primer_only_prefix']),
            'unresolved_reasons': ';'.join(sorted(reasons))}


def donor_fields():
    fields = ['read_id','mate','cleaned_raw_start_1based','cleaned_raw_end_1based',
              'junction_type','candidate_alignments','left_independent_matches',
              'right_independent_matches','unaligned_gap_bp','alignment_overlap_bp',
              'reference_rewind_bp']
    for side in ('left','right'):
        fields.extend(side+s for s in ('_read_start_1based','_read_end_1based',
            '_donor_start_1based','_donor_end_1based','_strand','_identity_percent',
            '_aligned_read','_aligned_donor'))
    return fields + ['join_subtype','accepted_host_interval_gap_overlap_bp']


HOST_FIELDS = ['read_id','mate','reference','position_1based','mapq','cigar','strand',
    'supplementary','accepted','identity_percent','independent_host_matches','match_method',
    'query_start_1based','query_end_1based','raw_query_start_1based','raw_query_end_1based','reasons']


def _metrics(pair_rows):
    counts = collections.Counter(r['classification'] for r in pair_rows)
    short = [r for r in pair_rows if r['classification']==CLASSES[2]]
    flags = ['both_supported','same_read_both_supported','host_mapping_supported',
             'donor_rearrangement_supported','abutting_or_overlapping_donor_anchors',
             'donor_gap_with_host_interval_candidate','donor_gap_other_candidate','R1_short_pattern']
    return {
        'input_pairs':len(pair_rows),
        'class_counts':{c:counts[c] for c in CLASSES},
        'both_supported_pairs':sum(int(r['both_supported']) for r in pair_rows),
        'host_supported_pairs':sum(int(r['host_mapping_supported']) for r in pair_rows),
        'donor_supported_pairs':sum(int(r['donor_rearrangement_supported']) for r in pair_rows),
        'short_pattern_pairs_any_class':sum(int(r['R1_short_pattern']) for r in pair_rows),
        'short_adapter_confirmed_pairs':sum(int(r['R1_adapter_confirmed']) for r in short),
        'short_polyA_only_pairs':sum(not int(r['R1_adapter_confirmed']) for r in short),
        'short_adapter_exact_pairs':sum(int(r['R1_adapter_exact']) for r in short),
        'short_primer_only_pairs':sum(int(r['R1_primer_only_prefix']) for r in short),
        'supplementary_flags':{f:sum(int(r[f]) for r in pair_rows) for f in flags},
    }


def write_chunk_outputs(records, raw_r1, donor, sample, side, outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    if set(records['R1']) != set(records['R2']) or not records['R1']:
        raise ValueError('Chunk requires nonempty, exactly matched R1/R2 IDs')
    pair_rows, read_rows, donor_rows, host_rows = [], [], [], []
    for rid, r1 in records['R1'].items():
        raw = raw_r1[rid]
        seq = raw[1] if not isinstance(raw, str) else raw
        pat = short_product_pattern(seq, donor, side)
        ends = [r1, records['R2'][rid]]
        pair_rows.append({'read_id':rid, **classify_pair(ends, pat)})
    for mate, group in records.items():
        for r in group.values():
            row = core.qc_row(r)
            row.update(donor_donor_supported=int(r['dd']), host_mapping_supported=int(r['host']),
                both_supported=int(r['dd'] and r['host']), donor_class=r['donor_class'],
                junction_type=r['donor_details'].get('junction_type',''),
                donor_masked_bases=len(r['donor_mask']),
                host_best_independent_matches=max((h.get('independent_host_matches',0) for h in r['host_hits']),default=0))
            read_rows.append(row)
            if r['dd']:
                subtype, overlap = donor_join_subtype(r)
                donor_rows.append({'read_id':r['read_id'],'mate':mate,
                    'cleaned_raw_start_1based':r['raw_start0']+1,'cleaned_raw_end_1based':r['raw_end0'],
                    **r['donor_details'],'join_subtype':subtype,
                    'accepted_host_interval_gap_overlap_bp':overlap})
            for h in r['host_hits']:
                host_rows.append({'read_id':r['read_id'],'mate':mate,**h,'reasons':';'.join(h['reasons'])})
    core.write_tsv(outdir/FILES[0],pair_rows,list(pair_rows[0]))
    core.write_tsv(outdir/FILES[1],read_rows,list(read_rows[0]))
    core.write_tsv(outdir/FILES[2],donor_rows,donor_fields())
    core.write_tsv(outdir/FILES[3],host_rows,HOST_FIELDS)
    metrics = _metrics(pair_rows)
    if sum(metrics['class_counts'].values()) != metrics['input_pairs']:
        raise ValueError('Four-class denominator mismatch')
    core.write_json(outdir/'chunk_metrics.json',metrics)
    return metrics


def _csv(path, rows, fields):
    tmp = path.with_name(path.name+'.tmp')
    with tmp.open('w',encoding='utf-8-sig',newline='') as f:
        w = csv.DictWriter(f,fieldnames=fields)
        w.writeheader(); w.writerows(rows)
    tmp.replace(path)


def _merge_tsv(paths, output):
    tmp = output.with_name(output.name+'.tmp')
    header = None
    n = 0
    with gzip.open(tmp,'wt',encoding='utf-8',newline='') as out:
        for path in paths:
            with path.open(encoding='utf-8',newline='') as src:
                current = src.readline()
                if not current:
                    raise ValueError('Empty TSV: '+str(path))
                if header is None:
                    header = current
                    out.write(header)
                elif current != header:
                    raise ValueError('Inconsistent TSV headers: '+str(path))
                for line in src:
                    if not line.strip():
                        raise ValueError('Empty TSV data row: '+str(path))
                    out.write(line); n += 1
    tmp.replace(output)
    return n


def aggregate_outputs(chunk_dirs, outdir, metadata):
    """Produce numeric source tables; no biological P-values or inferred n."""
    outdir = Path(outdir)
    n = int(metadata['total_pairs'])
    if n <= 0 or not chunk_dirs:
        raise ValueError('No pairs to aggregate')
    counts = collections.Counter()
    extras = collections.Counter()
    flags = collections.Counter()
    computed_n = 0
    for d in chunk_dirs:
        m = json.loads((Path(d)/'chunk_metrics.json').read_text())
        with (Path(d)/FILES[0]).open() as f:
            pairs = list(csv.DictReader(f,delimiter='\t'))
        actual = _metrics(pairs)
        if actual != m:
            raise ValueError('Chunk metrics differ from pair rows: '+str(d))
        if sum(m['class_counts'].values()) != m['input_pairs']:
            raise ValueError('Invalid chunk denominator: '+str(d))
        computed_n += m['input_pairs']; counts.update(m['class_counts'])
        flags.update(m['supplementary_flags'])
        extras.update({k:v for k,v in m.items() if k not in ('class_counts','supplementary_flags')})
    if computed_n != n or sum(counts.values()) != n:
        raise ValueError('Full input count differs from classified count')
    merged = {}
    for name in FILES:
        merged[name] = _merge_tsv([Path(d)/name for d in chunk_dirs], outdir/(name+'.gz'))
    if merged[FILES[0]] != n or merged[FILES[1]] != 2*n:
        raise ValueError('Merged pair/read count mismatch')
    base = {'sample':metadata['sample'],'library_side':metadata['side'],
            'condition':metadata.get('condition',''),
            'biological_replicate':metadata.get('biological_replicate','')}
    common = {'total_pairs':n,'denominator':'all_input_read_pairs',
              'analysis_version':metadata.get('analysis_version',VERSION)}
    rows = [{**base,'category':c,'category_label':LABELS[c],'count':counts[c],
             **common,'fraction':counts[c]/n,'percent':100*counts[c]/n} for c in CLASSES]
    _csv(outdir/'figure_counts.csv',rows,list(rows[0]))
    wide = {**base,**common}
    for c in CLASSES:
        wide[c+'_count'] = counts[c]
        wide[c+'_percent'] = 100*counts[c]/n
    wide['host_supported_including_donor_class_count'] = extras['host_supported_pairs']
    wide['both_supported_count'] = extras['both_supported_pairs']
    _csv(outdir/'figure_wide.csv',[wide],list(wide))
    supp = []
    descriptions = {
        'both_supported':'Donor and host support anywhere in the pair; included in donor class',
        'same_read_both_supported':'Both evidence types in at least one sequenced end',
        'host_mapping_supported':'All host-supported pairs including donor class',
        'donor_rearrangement_supported':'All donor-junction/rearrangement-supported pairs',
        'abutting_or_overlapping_donor_anchors':'At least one selected junction with no positive gap; not proof of integration',
        'donor_gap_with_host_interval_candidate':'Accepted host independent-match span lies within selected donor-anchor gap',
        'donor_gap_other_candidate':'Selected donor-anchor positive gap without contained accepted host span',
        'R1_short_pattern':'R1 short-product pattern irrespective of exclusive class',
    }
    for key, val in flags.items():
        supp.append({**base,'metric':key,'count':val,'denominator_count':n,
            'percent':100*val/n,'scope':'all_input_pairs_nonexclusive',
            'description':descriptions.get(key,key)})
    for key in ('short_adapter_confirmed_pairs','short_polyA_only_pairs','short_adapter_exact_pairs','short_primer_only_pairs'):
        den = counts[CLASSES[2]]
        supp.append({**base,'metric':key,'count':extras[key],'denominator_count':den,
            'percent':100*extras[key]/den if den else '',
            'scope':'short_donor_primer_candidate_class','description':key})
    _csv(outdir/'supplementary_counts.csv',supp,list(supp[0]))
    text = methods_text(metadata)
    (outdir/'methods.txt').write_text(text,encoding='utf-8')
    invariant = {'all_input_pairs_classified':computed_n==n,
                 'exclusive_classes_sum_to_total':sum(counts.values())==n,
                 'two_read_evidence_rows_per_pair':merged[FILES[1]]==2*n,
                 'host_support_reconciles':extras['host_supported_pairs']==counts[CLASSES[1]]+extras['both_supported_pairs'],
                 'short_confidence_subgroups_reconcile':extras['short_adapter_confirmed_pairs']+extras['short_polyA_only_pairs']==counts[CLASSES[2]]}
    if not all(invariant.values()):
        raise ValueError('Final count invariant failed: '+repr(invariant))
    return {'input_pairs':n,'total_pairs':n,'class_counts':{c:counts[c] for c in CLASSES},
            'class_percentages':{c:100*counts[c]/n for c in CLASSES},
            'both_supported_pairs':extras['both_supported_pairs'],
            'host_supported_pairs':extras['host_supported_pairs'],
            'short_adapter_confirmed_pairs':extras['short_adapter_confirmed_pairs'],
            'short_polyA_only_pairs':extras['short_polyA_only_pairs'],
            'validation':invariant,'evidence_row_counts':merged,
            'figure_files':['figure_counts.csv','figure_wide.csv','supplementary_counts.csv'],
            'interpretation':'Amplified sequencing-library read-pair composition; not cells, original molecules, independent insertions or proven concatemer integration.'}


def methods_text(m):
    return f'''RAISING four-class full-input analysis — method record v{VERSION}

Sample: {m['sample']}
Library direction: {m['side']}
Condition: {m.get('condition','')}
Biological replicate: {m.get('biological_replicate','')}
Number of input read pairs: {m['total_pairs']}
Denominator: every valid paired FASTQ record from all selected input files.
No random subsampling. No PCR deduplication or UMI collapse. Empty/short cleaned
reads remain in the denominator and can be unresolved. Malformed FASTQ, missing
mates, mismatched mate IDs, or duplicate IDs stop the run rather than being dropped.

Technical processing
The provided 800-bp donor and the 3-prime P3 primer are used. A leading R1 P3
or its reverse complement is recognized at a start offset of 0-12 nt with <=2
substitutions. R1 polyA (>=15) is trimmed only when followed by the complete
23-nt adapter TACTGCTTGACTGACCTGCTGTC with <=2 substitutions (including one-base
polyA boundary ambiguity). R2 leading reverse-complement adapter and adjacent
polyT (>=10) are trimmed. The 3-prime edge is trimmed using a 5-base mean Q20
window. PolyA alone is not trimmed. Details are preserved per read.

Donor junction/rearrangement evidence
Exact identical cleaned sequences within each chunk are searched once by
BLASTN (blastn-short, word_size 7, dust no, soft_masking false, evalue 100,
max_hsps 50). Results are returned to every original read without collapsing
their counts. Donor HSP identity >=95% and >=20 exact matches define the donor
query mask. A selected junction requires >=30 independently matching bases
on each side, excluding the P3 reference region and overlapping query bases;
identity >=95%, query gap <=50 bp and overlap <=20 bp. Same-direction pairs
must rewind >=20 bp on the reference. A continuous single-donor explanation
covering the same interval is rejected as junction evidence. Junction anchors
do not have an additional per-base Q20 criterion in this frozen screen.
Only the best qualifying junction per read is exported; this is not an
exhaustive catalog of all junctions. Positive evidence may reflect a concatemer,
partial rearrangement, circular donor, or PCR chimera, not necessarily integration.

Host evidence
Each cleaned end is independently aligned against the combined human + donor
reference using BWA-MEM -a -Y -T 20 -t <threads>. Read-pair insert-size constraints
are not used. Unmapped, secondary and donor-contig alignments are not accepted;
supplementary human alignments may contribute. Accepted evidence requires
MAPQ >=20 (255 excluded), conservative identity >=95%, and >=30 Q20 exact
matching query bases outside donor HSP intervals and homopolymers >=10 bases.
Identity = 100*(1-NM/(M+equals+X+I)); deletions count in NM. MD tags identify exact
matches, with a conservative NM fallback. Mapping supports a host-like sequence,
not proof of a donor-host chromosomal junction. Short 18-29-bp hits are not an
exhaustive hostcutoff18 reanalysis.

Short donor/primer product pattern (R1 only)
The prefix before the FIRST maximal A run of >=15 bases must be >=15 bases
and exactly match contiguous donor sequence starting at the P3 start or up to
2 bases downstream. It includes primer-only prefixes. The category is a
sequence-pattern candidate, not a call of residual donor, monomer, physical
donor end, or deletion length. Adapter-confirmed (<=2 substitutions at the
adjacent adapter) and polyA-only candidates are exported separately. Full exact
adapter matches are also counted. Uncertainty at an A boundary is not resolved.

Exclusive class priority
1 donor_rearrangement_supported: at least one end has donor junction evidence.
2 host_mappable_other: otherwise at least one end has accepted host evidence.
3 short_donor_primer_candidate: otherwise R1 has the above short-product pattern.
4 other_unresolved: all remaining pairs, not proven non-concatemer.
Both-support flags, donor junction orientation, and short-pattern confidence
are retained separately. A host independent-match coordinate span contained
within a selected positive donor-anchor gap is a supplementary D-H-D candidate
flag, not independent verification of a copied host fragment.

Outputs and interpretation
figure_counts.csv contains one row per exclusive category, with zero categories
retained, integer counts, denominator, fraction (0-1), and percent (0-100).
figure_wide.csv contains one sample row. supplementary_counts.csv contains
nonexclusive or explicitly subgrouped counts; do not sum its rows as a four-class
composition. Evidence TSVs and chunk results retain read-level auditability.
run_summary.json must have status=complete and all validation checks true before
using final CSVs. Results quantify the amplified sequencing library. Statistical
comparisons must use independent experimental replicates, not read pairs as n.
5-prime and 3-prime libraries from the same sample are not separate biological
replicates. This release's executable supports 3-prime libraries only.
MH lengths, terminal deletion lengths, independent insertion-site counts and
PCR-chimera exclusion have not been determined by this classification run.

Software versions (exact runtime record)
{json.dumps(m.get('tool_versions',{}),ensure_ascii=False,indent=2)}

Configuration
{json.dumps(m.get('config',{}),ensure_ascii=False,indent=2,default=str)}

Input identity/checksums and full runtime provenance are in run_summary.json.
'''
