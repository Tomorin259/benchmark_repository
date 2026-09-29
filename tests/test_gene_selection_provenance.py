import ast
import csv
import hashlib
import json
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SELECTION = ROOT / 'evaluation/gene/selection'


class GeneSelectionProvenanceTests(unittest.TestCase):
    def test_preserved_sources_match_checksums_and_rows(self):
        manifest = json.loads((SELECTION / 'provenance.json').read_text())
        for record in manifest['historical_files']:
            path = SELECTION / record['snapshot']
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), record['sha256'])
            self.assertTrue(record['source_id'].startswith('artifact://velocity/'))
            if path.suffix == '.csv':
                with path.open() as stream:
                    count = sum(1 for _ in csv.reader(stream))
                self.assertEqual(count - int(record['csv_has_header']), record['csv_data_rows'])
        upstream = manifest['upstream']
        self.assertEqual(hashlib.sha256((SELECTION / upstream['local_entry']).read_bytes()).hexdigest(), upstream['sha256'])
        self.assertTrue((SELECTION / 'upstream/LICENCE.txt').is_file())

    def test_frozen_selection_matches_evaluation_and_candidate_evidence(self):
        tree = ast.parse((ROOT / 'evaluation/gene/run_reference30_gene_cbdir.py').read_text())
        configs = next(ast.literal_eval(n.value) for n in tree.body
                       if isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'REAL_CONFIGS' for t in n.targets))
        with (SELECTION / 'selected_genes.tsv').open() as stream:
            rows = list(csv.DictReader(stream, delimiter='\t'))
        pairs = [(r['dataset'], r['gene']) for r in rows]
        self.assertEqual(len(pairs), len(set(pairs)))
        self.assertEqual(set(pairs), {(d, gene) for d in {r['dataset'] for r in rows} for gene in configs[d]['genes']})
        for row in rows:
            for candidate in row['candidate_evidence'].split(';'):
                with (SELECTION / candidate).open() as stream:
                    data = list(csv.DictReader(stream))
                key = 'gene' if 'gene' in data[0] else next(iter(data[0]))
                self.assertIn(row['gene'], {r[key] for r in data})
            self.assertEqual(row['manual_reason_status'], 'not_recorded')

    def test_reconstruction_entry_points_compile_and_have_portable_paths(self):
        for name in ['reproduce_murk.py', 'reproduce_branching.py', 'compare_rankings.py']:
            source = (SELECTION / name).read_text()
            compile(source, name, 'exec')
            self.assertNotIn('/mnt/', source)
            self.assertNotIn('/Users/', source)

    def test_reconstruction_summary_matches_preserved_comparisons(self):
        summary = json.loads((SELECTION / 'reconstruction/summary.json').read_text())
        self.assertEqual(summary['frozen_evaluation_dataset_gene_pairs'], 59)
        self.assertEqual(len(summary['attempts']), 3)
        self.assertEqual(summary['final_selection']['per_gene_decision_logs'], 'not_retained')
        for attempt in summary['attempts']:
            self.assertTrue((SELECTION / attempt['manifest']).is_file())
            for check in attempt['checks']:
                comparison = json.loads((SELECTION / check['comparison']).read_text())
                for field in ['n_historical', 'n_rerun', 'n_common', 'row_order_exact',
                              'tie_aware_primary_score_ranks_equal', 'top50_overlap']:
                    self.assertEqual(check[field], comparison[field])
                for role in ['historical', 'rerun']:
                    path = SELECTION / check[role]
                    self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(),
                                     comparison[role + '_sha256'])
                if check['status'] == 'exact_values_and_tie_aware_ranking':
                    self.assertTrue(check['candidate_set_equal'])
                    self.assertTrue(check['common_gene_numeric_fields_exact'])
                    self.assertTrue(check['tie_aware_primary_score_ranks_equal'])

    def test_recovered_dentategyrus_groups_and_ranks_match_snapshots(self):
        base = ast.parse((SELECTION / 'historical/scripts/branching_DentateGyrus_lamanno_20251216.py').read_text())
        pro = ast.parse((SELECTION / 'historical/scripts/branching_DentateGyrus_lamanno_pro_20260107.py').read_text())
        entry = ast.parse((SELECTION / 'reproduce_branching.py').read_text())
        def assigned(tree, name):
            return next(ast.literal_eval(n.value) for n in ast.walk(tree)
                        if isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id == name for t in n.targets))
        granule = next(ast.literal_eval(n.value)[0] for n in ast.walk(entry)
                       if isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'tasks' for t in n.targets)
                       and ast.literal_eval(n.value)[0][0] == 'Granule_vs_Pyramidal')
        neuron = next(ast.literal_eval(n.args[0]) for n in ast.walk(entry)
                      if isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute)
                      and n.func.attr == 'append' and n.args and isinstance(n.args[0], ast.Tuple)
                      and ast.literal_eval(n.args[0])[0] == 'Neuron_vs_Glia')
        self.assertEqual(neuron[1], assigned(base, 'cloud_neuron'))
        self.assertEqual(neuron[2], assigned(base, 'cloud_glia'))
        self.assertEqual(neuron[3:], (None, None))
        self.assertEqual(granule[1], assigned(base, 'cloud_granule'))
        self.assertEqual(granule[2], assigned(base, 'cloud_pyramidal'))
        self.assertEqual(granule[3], assigned(pro, 'rank_dict_granule'))
        self.assertEqual(granule[4], assigned(pro, 'rank_dict_pyramidal'))


if __name__ == '__main__':
    unittest.main()
