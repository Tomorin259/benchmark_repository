"""Exploratory case exports stay synchronized with the current score inputs."""

import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "supplement/gene-kinetic/analysis"


class GeneCaseResultTests(unittest.TestCase):
    def test_committed_analysis_replays_from_current_inputs(self):
        with tempfile.TemporaryDirectory() as folder:
            subprocess.run(
                [sys.executable, str(ANALYSIS / "gene_kinetic_analysis.py"),
                 "--repo-root", str(ROOT), "--output-dir", folder],
                check=True, capture_output=True, text=True,
            )
            outputs = list(Path(folder).iterdir())
            self.assertGreater(len(outputs), 10)
            for path in outputs:
                retained = ANALYSIS / path.name
                if path.name == "provenance.json":
                    self.assertEqual(json.loads(path.read_text())["inputs_sha256"],
                                     json.loads(retained.read_text())["inputs_sha256"])
                elif path.suffix == ".csv":
                    with path.open() as actual_file, retained.open() as expected_file:
                        actual = list(csv.reader(actual_file))
                        expected = list(csv.reader(expected_file))
                    self.assertEqual(len(actual), len(expected), path.name)
                    self.assertEqual(actual[0], expected[0], path.name)
                    for actual_row, expected_row in zip(actual[1:], expected[1:]):
                        self.assertEqual(len(actual_row), len(expected_row))
                        for value, reference in zip(actual_row, expected_row):
                            if value != reference:
                                try:
                                    difference = abs(float(value) - float(reference))
                                except ValueError:
                                    self.assertEqual(value, reference, path.name)
                                else:
                                    self.assertLessEqual(difference, 1e-12, path.name)
                else:
                    self.assertEqual(path.read_bytes(), retained.read_bytes(), path.name)

    def test_portrait_annotations_match_current_gene_scores(self):
        spec = importlib.util.spec_from_file_location("case_analysis", ANALYSIS / "gene_kinetic_analysis.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        manifests = module.load_manifests(str(ROOT))
        _, cells, all_cells = module.load_cell_tables(str(ROOT))
        mapping = module.build_source_mapping(manifests, module.load_audits(str(ROOT)))
        denominators, targets = module.make_branch_denominators(all_cells)
        rows = module.make_observations(manifests, module.load_gene_records(str(ROOT)),
                                        cells, mapping, denominators, targets)
        self.assertEqual(len(rows), 252)
        for row in rows:
            with self.subTest(case=row["case"], method=row["method"], edge=(row["source"], row["target"])):
                if row["gene_cbdir"] is None:
                    self.assertIsNone(row["manifest_edge_score"])
                else:
                    self.assertIsNotNone(row["manifest_edge_score"])
                    self.assertLessEqual(abs(row["manifest_edge_score"] - row["gene_cbdir"]), 0.0005000001)

    def test_analysis_sources_and_outputs_match_provenance(self):
        provenance = json.loads((ANALYSIS / "provenance.json").read_text())
        for key, base in [("inputs_sha256", ROOT), ("outputs_sha256", ANALYSIS)]:
            for relative, expected in provenance[key].items():
                self.assertEqual(hashlib.sha256((base / relative).read_bytes()).hexdigest(), expected, relative)


if __name__ == "__main__":
    unittest.main()
