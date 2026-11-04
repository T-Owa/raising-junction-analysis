"""Synthetic tests for counting, metadata, missing values, and input validation.

No experimental reads or sample identifiers are included in these fixtures.
"""
import contextlib
import csv
import gzip
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / 'tools' / 'prepare_figure_data.py'
SPEC = importlib.util.spec_from_file_location('public_figure_data', SCRIPT)
MOD = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MOD)


def write_table(path, fields, rows):
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def fixture(directory, include_candidates=True):
    """Construct two arbitrary named libraries, one with no junction candidates."""
    directory.mkdir()
    sample = 'experiment-A.left'
    groups, annotations, comparisons = [], [], []
    configs = [
        # id, support, simple, prior-DD, complex, MH, loss-min, insert, compartment
        ('J1', 3, 3, 0, 0, 0, 0, '', 'nuclear_or_unplaced_reference'),
        ('J2', 2, 1, 1, 0, 2, 1, '', 'nuclear_or_unplaced_reference'),
        ('J3', 2, 0, 0, 2, 0, 0, 'GC', 'nuclear_or_unplaced_reference'),
        ('J4', 4, 4, 0, 0, 1, 0, '', 'mitochondrial_reference'),
        ('J5', 0, 0, 0, 0, 0, 2, '', 'nuclear_or_unplaced_reference'),
        ('J6', 1, 1, 0, 0, 1, 0, '', 'nuclear_or_unplaced_reference'),
    ] if include_candidates else []
    for index, (event, support, simple, dd, complex_count, mh, loss, insertion, compartment) in enumerate(configs):
        meta = dict(sample=sample, library_side='3p', condition='', biological_replicate='')
        group_id = 'E' + event[1:]
        review = 1 if event == 'J5' else 0
        groups.append({**meta, 'equivalent_template_group': group_id,
                       'representative_template_id': event, 'template_ids': event,
                       'reference_compartment': compartment, 'supported_pairs': support,
                       'review_pairs': review, 'without_prior_dd_or_complex_host_pairs': simple,
                       'prior_dd_without_complex_host_pairs': dd, 'complex_host_pairs': complex_count})
        row = dict.fromkeys(MOD.ANNOTATION_KEEP, '')
        row.update(dict.fromkeys(MOD.ANNOTATION_NUMERIC, 0))
        row.update({**meta, 'event_id': event, 'equivalent_template_group': group_id,
                    'reference_compartment': compartment, 'chrom': 'MT' if event == 'J4' else '1',
                    'host_strand': '+', 'host_boundary0': 100 + index * 10,
                    'host_boundary_min0': 100 + index * 10,
                    'host_boundary_max0': 100 + index * 10 + mh,
                    'donor_deletion_min_bp': loss, 'donor_deletion_max_bp': loss + mh,
                    'mh_length_bp': mh, 'mh_sequence': 'A' * mh,
                    'insertion_length_bp': len(insertion), 'insertion_sequence': insertion,
                    'pair_support': support, 'simple_pair_support': simple,
                    'prior_dd_pair_support': dd, 'complex_host_pair_support': complex_count,
                    'is_group_representative': 1,
                    'supported_pair_count_for_identical_sequence_group': support,
                    'supported_pair_count_once_per_group': support,
                    'review_pair_count_for_identical_sequence_group': review,
                    'template_reference_mh_length_bp': mh,
                    'template_reference_donor_deletion_min_bp': loss,
                    'template_reference_donor_deletion_max_bp': loss + mh,
                    'mh_for_analysis_bp': '' if insertion else mh,
                    'mh_assessment': ('not_assessed_insertion_interfaces' if insertion
                                      else 'direct_junction_reference_annotation')})
        annotations.append(row)
        comparison = dict.fromkeys(MOD.COMPARISON_KEEP, '')
        comparison.update(sample=sample, template_id=event, reference_template_id=event,
                          relationship='strongest_template',
                          has_junction_proximal_difference='False',
                          differences_only_donor_internal_substitutions='False',
                          alignment_ambiguous='False', alignment_edit_bases=0,
                          substitutions=0, indel_bases=0, gap_events=0, operations_json='[]')
        comparisons.append(comparison)
    group_fields = [*dict(sample='', library_side='', condition='', biological_replicate=''),
                    'equivalent_template_group','representative_template_id','template_ids',
                    'reference_compartment',*MOD.GROUP_NUMERIC]
    write_table(directory / MOD.TABLE_NAMES[0], group_fields, groups)
    write_table(directory / MOD.TABLE_NAMES[1], MOD.ANNOTATION_KEEP, annotations)
    write_table(directory / MOD.TABLE_NAMES[2], MOD.COMPARISON_KEEP, comparisons)
    libraries = []
    for name, side in [(sample, '3p'), ('empty.right', '5p')]:
        have = name == sample and include_candidates
        libraries.append({
            'sample': name, 'library_side': side, 'status': 'complete',
            'original_candidate_patterns': 6 if have else 0,
            'identical_sequence_groups': 6 if have else 0,
            'supported_identical_sequence_groups': 5 if have else 0,
            'total_input_pairs': 30 if have else 10,
            'host_supported_pairs': 15 if have else 0,
            **dict.fromkeys(MOD.STATUS_KEYS, 0),
            'template_supported': 12 if have else 0,
            'junction_variant_review': 1 if have else 0,
            'no_nearby_template': 2 if have else 0,
            'supported_pairs_by_reference_compartment': {
                'nuclear_or_unplaced_reference': 8 if have else 0,
                'mitochondrial_reference': 4 if have else 0},
        })
    summary = {
        'status': 'complete', 'completed_library_count': 2, 'selected_library_count': 2,
        'libraries': libraries,
        'merged_outputs': {name + '.gz': {'rows': len(groups)} for name in MOD.TABLE_NAMES},
    }
    (directory / 'batch_summary.json').write_text(json.dumps(summary))
    return summary


class FigureDataTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.source = self.root / 'input'
        self.output = self.root / 'output'
        self.summary = fixture(self.source)

    def tearDown(self):
        self.temp.cleanup()

    def run_export(self, *extra):
        with contextlib.redirect_stdout(io.StringIO()):
            return MOD.main(['--input-dir',str(self.source),'--output-dir',str(self.output),*extra])

    def output_rows(self, name):
        with (self.output / name).open(encoding='utf-8-sig') as handle:
            return list(csv.DictReader(handle))

    def mutate_table(self, name, function):
        path = self.source / name
        with path.open() as handle:
            reader=csv.DictReader(handle,delimiter='\t'); records=list(reader); fields=reader.fieldnames
        function(records)
        write_table(path,fields,records)

    def test_denominators_support_subsets_and_na(self):
        self.run_export()
        data = {(r['sample'],r['subset'],int(r['min_support'])):r
                for r in self.output_rows('figure_values.csv')}
        primary = data['experiment-A.left','nuclear_all',1]
        self.assertEqual(primary['N_patterns'],'4')
        self.assertEqual(primary['N_direct'],'3')
        self.assertAlmostEqual(float(primary['MH_ge2_fraction']),1/3)
        self.assertAlmostEqual(float(primary['donor_loss_min_gt0_fraction']),1/3)
        self.assertEqual(float(primary['insertion_candidate_fraction']),0.25)
        self.assertEqual(data['experiment-A.left','nuclear_with_simple',2]['N_patterns'],'1')
        self.assertEqual(data['experiment-A.left','nuclear_only_simple',1]['N_patterns'],'2')
        rows = {r['event_id']:r for r in self.output_rows('candidate_annotations.csv')}
        self.assertEqual(rows['J3']['mh_for_analysis_bp'],'')
        self.assertEqual(primary['condition'],'')
        self.assertEqual(primary['series_code'],'')

    def test_zero_candidate_library_has_zero_counts_and_missing_fractions(self):
        self.run_export()
        zero = [r for r in self.output_rows('figure_values.csv') if r['sample']=='empty.right']
        self.assertEqual(len(zero),6)
        self.assertTrue(all(r['N_patterns']=='0' and r['MH0_fraction']==''
                            and r['insertion_candidate_fraction']=='' for r in zero))

    def test_entire_batch_without_candidates_has_output_headers(self):
        empty=self.root/'empty_input';fixture(empty,include_candidates=False);self.source=empty
        result=self.run_export()
        self.assertEqual(result['totals']['candidate_groups'],0)
        self.assertEqual(len(self.output_rows('figure_values.csv')),12)
        self.assertEqual(self.output_rows('candidate_annotations.csv'),[])
        self.assertEqual(self.output_rows('alternative_review.csv'),[])
        self.assertIn('event_id',(self.output/'candidate_annotations.csv').read_text())

    def test_compressed_tables_match_uncompressed(self):
        self.run_export()
        previous=(self.output/'figure_values.csv').read_bytes()
        for name in MOD.TABLE_NAMES:
            path=self.source/name
            with gzip.open(str(path)+'.gz','wb') as handle:handle.write(path.read_bytes())
            path.unlink()
            self.summary['merged_outputs'][name+'.gz']['sha256']=MOD.digest(Path(str(path)+'.gz'))
        (self.source/'batch_summary.json').write_text(json.dumps(self.summary))
        self.run_export()
        self.assertEqual((self.output/'figure_values.csv').read_bytes(),previous)

    def test_duplicate_compressed_and_plain_inputs_fail(self):
        path=self.source/MOD.TABLE_NAMES[0]
        with gzip.open(str(path)+'.gz','wb') as handle:handle.write(path.read_bytes())
        with self.assertRaisesRegex(ValueError,'Expected exactly one'):
            self.run_export()

    def test_explicit_metadata_fills_blanks_without_parsing_names(self):
        path=self.root/'metadata.csv'
        path.write_text('sample,library_side,condition,biological_replicate\nexperiment-A.left,3p,treatment-A,replicate-X\n')
        self.run_export('--metadata',str(path))
        row=self.output_rows('figure_values.csv')[6]  # sorted order: empty.right first
        self.assertEqual(row['condition'],'treatment-A')
        self.assertEqual(row['biological_replicate'],'replicate-X')
        self.assertEqual(row['replicate_label_from_name'],'')

    def test_existing_table_metadata_is_retained(self):
        self.mutate_table(MOD.TABLE_NAMES[0],lambda rows:[r.update(condition='control',biological_replicate='B1') for r in rows])
        self.run_export()
        self.assertTrue(all(r['condition']=='control' for r in self.output_rows('candidate_annotations.csv')))

    def test_conflicting_metadata_fails(self):
        self.mutate_table(MOD.TABLE_NAMES[0],lambda rows:rows[0].update(condition='control'))
        self.mutate_table(MOD.TABLE_NAMES[1],lambda rows:rows[0].update(condition='treated'))
        with self.assertRaisesRegex(ValueError,'Conflicting condition'):
            self.run_export()

    def test_support_subtotals_must_reconcile(self):
        self.mutate_table(MOD.TABLE_NAMES[0],lambda rows:rows[0].update(supported_pairs='9'))
        with self.assertRaises(ValueError):self.run_export()

    def test_multiple_member_groups_are_explicitly_rejected(self):
        self.mutate_table(MOD.TABLE_NAMES[0],lambda rows:rows[0].update(template_ids='J1;J2'))
        with self.assertRaisesRegex(ValueError,'single-member groups'):
            self.run_export()

    def test_insertion_mh_must_be_unassessed(self):
        self.mutate_table(MOD.TABLE_NAMES[1],lambda rows:rows[2].update(mh_for_analysis_bp='0'))
        with self.assertRaises(ValueError):self.run_export()

    def test_recorded_batch_rows_must_match(self):
        self.summary['merged_outputs'][MOD.TABLE_NAMES[0]+'.gz']['rows']=7
        (self.source/'batch_summary.json').write_text(json.dumps(self.summary))
        with self.assertRaisesRegex(ValueError,'Merged row count mismatch'):
            self.run_export()


if __name__ == '__main__':
    unittest.main()
