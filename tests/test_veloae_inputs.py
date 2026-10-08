"""Frozen PCA source coverage and projection identities for Supplementary Fig. 11."""
import csv
import hashlib
import json
from pathlib import Path
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
MODULE = ROOT / 'supplement/veloae'


class VeloAEInputTests(unittest.TestCase):
    def test_paired_pca_table_matches_frozen_hash_and_paper_means(self):
        path = MODULE / 'inputs/pca_comparison.csv'
        retained = json.loads((MODULE / 'results/pca_summary_manifest.json').read_text())
        self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), retained['input_sha256'])
        with path.open() as handle:
            rows = list(csv.DictReader(handle))
        keys = [(r['dataset'], r['method'], r['matrix_key']) for r in rows]
        self.assertEqual(len(rows), 84)
        self.assertEqual(len(keys), len(set(keys)))
        for method in ['veloae', 'scvelo_stochastic', 'TIvelo']:
            for matrix in ['Ms', 'velocity']:
                subset = [r for r in rows if r['method'] == method and r['matrix_key'] == matrix]
                self.assertEqual(len({r['dataset'] for r in subset}), 14)
        for matrix, expected in [('Ms', (95.7, 70.0)), ('velocity', (84.8, 52.7))]:
            subset = [r for r in rows if r['method'] == 'veloae' and r['matrix_key'] == matrix]
            means = tuple(round(sum(float(r[k]) for r in subset) / 14 * 100, 1)
                          for k in ['raw_pc1', 'standardized_pc1'])
            self.assertEqual(means, expected)

    def test_all_projection_inputs_match_retained_geometry_manifest(self):
        sources = json.loads((MODULE / 'inputs/manifest.json').read_text())['files']
        expected = {p['npz_sha256']: p for p in json.loads((MODULE / 'results/manifest.json').read_text())['panels']}
        projections = [r for r in sources if r['path'].endswith('.npz')]
        self.assertEqual(len(projections), 14)
        for record in sources:
            path = MODULE / record['path']
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), record['sha256'])
            self.assertTrue(record['source'].startswith('artifact://velocity/'))
            if path.suffix != '.npz':
                continue
            with np.load(path, allow_pickle=True) as data:
                retained = expected[record['sha256']]
                self.assertEqual(data['X2'].shape, (retained['n_cells'], 2))
                self.assertEqual(data['V2'].shape, data['X2'].shape)
                self.assertEqual(len(data['labels']), retained['n_cells'])
                self.assertEqual(data['label_categories'].astype(str).tolist(), retained['categories'])
                np.testing.assert_allclose(data['explained_variance_ratio'], retained['explained_variance_ratio'])


if __name__ == '__main__':
    unittest.main()
