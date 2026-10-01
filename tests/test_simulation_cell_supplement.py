import csv
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
MODULE = REPO / 'evaluation/cell/simulation'
SCRIPT = MODULE / 'build_score_table.py'


class SimulationCellSupplementTests(unittest.TestCase):
    def test_frozen_exports_match_full_precision_results(self):
        with tempfile.TemporaryDirectory() as directory:
            subprocess.run([sys.executable, str(SCRIPT), '--check'],
                           cwd=directory, check=True, capture_output=True)
        with (MODULE / 'results/simulation_cell_scores_4dp.csv').open(newline='') as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 27)
        for row in rows:
            for metric in ('Consistency Score', 'CBDir', 'ICCoh', 'TransCosine', 'BATC'):
                self.assertRegex(row[metric], r'^-?\d+\.\d{4}$')

    def test_reexport_preserves_published_csv_and_latex(self):
        with tempfile.TemporaryDirectory() as directory:
            subprocess.run([sys.executable, str(SCRIPT), '--output-root', directory],
                           check=True, capture_output=True)
            for filename in ('simulation_cell_scores_4dp.csv', 'simulation_cell_scores.tex',
                             'table_manifest.json'):
                self.assertEqual((Path(directory) / filename).read_bytes(),
                                 (MODULE / 'results' / filename).read_bytes())

    def test_check_rejects_changed_csv_and_latex(self):
        with tempfile.TemporaryDirectory() as directory:
            for filename in ('simulation_cell_scores_4dp.csv', 'simulation_cell_scores.tex'):
                subprocess.run([sys.executable, str(SCRIPT), '--output-root', directory],
                               check=True, capture_output=True)
                path = Path(directory) / filename
                path.write_text(path.read_text().replace('0.9667', '0.9666', 1))
                result = subprocess.run(
                    [sys.executable, str(SCRIPT), '--output-root', directory, '--check'],
                    capture_output=True, text=True,
                )
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('Stale supplementary table export', result.stderr)

    def test_final_pdf_matches_recorded_export(self):
        manifest = json.loads((REPO / 'figures/suppfig12_composite_manifest.json').read_text())
        pdf = (REPO / manifest['composite_path']).read_bytes()
        self.assertTrue(pdf.startswith(b'%PDF-'))
        self.assertEqual(hashlib.sha256(pdf).hexdigest(), manifest['composite_sha256'])
        source = REPO / manifest['metric_source_path']
        self.assertEqual(hashlib.sha256(source.read_bytes()).hexdigest(),
                         manifest['metric_source_sha256'])


if __name__ == '__main__':
    unittest.main()
