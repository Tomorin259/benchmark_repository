"""Simulation input layers and named ODE reference arrays."""
import csv
import hashlib
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd

from evaluation.gene import config
from evaluation.gene.simulation.reference import materialize_clean_reference as make_reference


ROOT = Path(__file__).resolve().parents[1]


class CleanInputReleaseTests(unittest.TestCase):
    def test_input_and_source_identities_are_distinct(self):
        with (ROOT / 'data/input_release_manifest.tsv').open() as handle:
            inputs = [r for r in csv.DictReader(handle, delimiter='\t') if r['role'] == 'fixed_simulation']
        with (ROOT / 'data/source_provenance.tsv').open() as handle:
            sources = {r['current_artifact_id']: r for r in csv.DictReader(handle, delimiter='\t')
                       if r['role'] == 'simulation_source'}
        self.assertEqual(len(inputs), 3)
        for row in inputs:
            self.assertNotEqual(row['sha256'], sources[row['artifact_id']]['sha256'])
            self.assertTrue({'clean_spliced', 'clean_unspliced', 'velocity_s_true', 'velocity_u_true'}
                            .issubset(row['layers'].split(';')))
            self.assertNotEqual(row['public_url'], sources[row['artifact_id']]['url'])

    def test_clean_layers_materialize_without_a_simulator_and_unknown_bytes_fail(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            repository, artifacts = root / 'repo', root / 'artifacts'
            (repository / 'data').mkdir(parents=True)
            artifacts.mkdir()
            cases, rows = {}, []
            for scenario in ['murk', 'branching', 'all']:
                obj = ad.AnnData(np.zeros((2, 1), dtype=np.float32),
                                 obs=pd.DataFrame(index=['0', '1']), var=pd.DataFrame(index=['g']))
                for key, value in {'alpha_true': 5., 'beta_true': 1., 'gamma_true': 2.,
                                   'clean_spliced': 3., 'clean_unspliced': 1.,
                                   'velocity_s_true': -5., 'velocity_u_true': 4.}.items():
                    obj.layers[key] = np.full((2, 1), value)
                name = scenario + '.h5ad'
                obj.write_h5ad(artifacts / name)
                digest = hashlib.sha256((artifacts / name).read_bytes()).hexdigest()
                cases[scenario] = {'reference_source': name, 'genes': ['g']}
                rows.append({'path': name, 'sha256': digest})
            with (repository / 'data/input_file_manifest.tsv').open('w', newline='') as handle:
                writer = csv.DictWriter(handle, ['path', 'sha256'], delimiter='\t')
                writer.writeheader(); writer.writerows(rows)
            with patch.object(config, 'SIM_CONFIGS', cases):
                records = make_reference(repository, artifacts, root / 'reference')
                self.assertEqual(len(records), 3)
                for record in records:
                    with np.load(root / 'reference' / record['output'], allow_pickle=False) as saved:
                        np.testing.assert_array_equal(saved['velocity_s_true'], [[-5.], [-5.]])
                        np.testing.assert_array_equal(saved['velocity_u_true'], [[4.], [4.]])
                with (repository / 'data/input_file_manifest.tsv').open('w', newline='') as handle:
                    writer = csv.DictWriter(handle, ['path', 'sha256'], delimiter='\t')
                    writer.writeheader(); writer.writerows([{**r, 'sha256': '1' * 64} for r in rows])
                with self.assertRaisesRegex(ValueError, 'hash differs'):
                    make_reference(repository, artifacts, root / 'unknown')


if __name__ == '__main__':
    unittest.main()
