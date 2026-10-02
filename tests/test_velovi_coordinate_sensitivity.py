"""Cheap checks of Table S6 provenance and the coordinate transformation."""

import csv
import hashlib
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
MODULE = REPO / "supplement/velovi"
SCRIPT = MODULE / "run_coordinate_sensitivity.py"
spec = importlib.util.spec_from_file_location("velovi_diagnostic", SCRIPT)
diagnostic = importlib.util.module_from_spec(spec)
spec.loader.exec_module(diagnostic)


class VeloVICoordinateTests(unittest.TestCase):
    def test_source_contract_and_primary_benchmark(self):
        manifest = json.loads((MODULE / "results/manifest.json").read_text())
        rows = list(csv.DictReader((MODULE / "results/full_dataset_three_condition.csv").read_text().splitlines()))
        self.assertEqual(len(rows), 42)
        self.assertEqual(len({row["dataset"] for row in rows}), 14)
        primary = {
            row["dataset"]: float(row["transcosine"])
            for row in csv.DictReader((REPO / "evaluation/cell/results/cell_metrics.csv").read_text().splitlines())
            if row["method"] == "veloVI"
        }
        for item in manifest["inputs"]:
            selected = [row for row in rows if row["dataset"] == item["dataset"]]
            self.assertEqual({row["condition"] for row in selected}, set(diagnostic.CONDITIONS))
            for row in selected:
                self.assertEqual(row["status"], "ok")
                self.assertEqual(int(row["n_edges"]), len(item["cluster_edges"]))
                self.assertEqual(int(row["n_neighbors"]), 30)
                for key in ("h5ad_sha256", "moments_sha256", "neighbors_sha256"):
                    self.assertEqual(row[key], item[key])
            values = {row["condition"]: float(row["transcosine"]) for row in selected}
            self.assertGreater(values[diagnostic.CONDITIONS[2]], values[diagnostic.CONDITIONS[1]])
            self.assertAlmostEqual(values[diagnostic.CONDITIONS[1]], primary[item["dataset"]], delta=2e-8)
        for name, expected in manifest["outputs"].items():
            self.assertEqual(hashlib.sha256((MODULE / "results" / name).read_bytes()).hexdigest(), expected)
        self.assertEqual(hashlib.sha256(SCRIPT.read_bytes()).hexdigest(), manifest["script_sha256"])

    def test_table_regeneration(self):
        with tempfile.TemporaryDirectory() as directory:
            subprocess.run([sys.executable, str(SCRIPT), "--table-only", "--output-root", directory],
                           cwd=directory, check=True, capture_output=True)
            new = list(csv.DictReader((Path(directory) / "table_s6.csv").read_text().splitlines()))
            old = list(csv.DictReader((MODULE / "results/table_s6.csv").read_text().splitlines()))
            self.assertEqual(len(new), 15)
            for generated, stored in zip(new, old):
                self.assertEqual(generated["dataset"], stored["dataset"])
                for key in diagnostic.TABLE_FIELDS[1:]:
                    self.assertAlmostEqual(float(generated[key]), float(stored[key]), delta=1e-15)
            self.assertEqual([f"{float(new[-1][key]):.4f}" for key in diagnostic.TABLE_FIELDS[1:]],
                             ["0.1147", "0.0778", "0.2193", "0.1415"])
        rows = list(csv.DictReader((MODULE / "results/full_dataset_three_condition.csv").read_text().splitlines()))
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(ValueError):
                diagnostic.write_table_s6(rows + [rows[0]], Path(directory) / "duplicate.csv")
            with self.assertRaises(ValueError):
                diagnostic.write_table_s6(rows[:-1], Path(directory) / "missing.csv")

    def test_identifier_alignment_and_affine_derivative(self):
        try:
            import numpy as np
        except ImportError:
            self.skipTest("NumPy is required for the synthetic matrix test")
        native = np.asarray([[0., 1.], [1., 3.], [3., 5.]])
        expected = native * np.asarray([100., .1]) + np.asarray([7., -2.])
        shuffled = expected[[2, 0, 1]]
        aligned = diagnostic.align_reference_ms(shuffled, ["c", "a", "b"], ["a", "b", "c"])
        np.testing.assert_array_equal(aligned, expected)
        slope, intercept, r2 = diagnostic.fit_affine(native, aligned)
        np.testing.assert_allclose(slope, [100., .1])
        np.testing.assert_allclose(intercept, [7., -2.])
        np.testing.assert_allclose(r2, 1.)
        # The derivative uses slopes only; translations must not enter velocity.
        np.testing.assert_allclose(np.asarray([[1., 2.]]) * slope, [[100., .2]])
        with self.assertRaises(ValueError):
            diagnostic.align_reference_ms(expected, ["a", "a", "c"], ["a", "b", "c"])
        with self.assertRaises(ValueError):
            diagnostic.align_reference_ms(expected, ["a", "b", "c"], ["a", "b", "d"])

    def test_online_table_s6_four_decimal_values(self):
        published = {
            "pancreas": ("0.1131", "0.0700", "0.2927", "0.2227"),
            "gastrulation": ("0.0511", "0.0603", "0.3473", "0.2870"),
            "forebrain": ("0.0904", "0.0661", "0.1528", "0.0867"),
            "DentateGyrus_lamanno": ("0.0588", "0.0537", "0.1709", "0.1172"),
            "dentategyrus": ("0.0955", "0.0626", "0.1108", "0.0482"),
            "retina": ("0.1592", "0.1418", "0.3040", "0.1622"),
            "mousebrain": ("0.1316", "0.1077", "0.3415", "0.2338"),
            "hematopoiesis": ("0.1271", "0.0683", "0.1455", "0.0772"),
            "bonemarrow": ("0.0352", "0.0600", "0.3338", "0.2738"),
            "rpe1": ("0.1117", "0.0860", "0.1252", "0.0393"),
            "organoid": ("0.0852", "0.0386", "0.1172", "0.0785"),
            "neuron": ("0.1216", "0.1059", "0.2461", "0.1402"),
            "capofibro": ("0.2541", "0.0688", "0.2433", "0.1745"),
            "rpe1_gen": ("0.1719", "0.1000", "0.1398", "0.0399"),
            "Macro-average": ("0.1147", "0.0778", "0.2193", "0.1415"),
        }
        rows = csv.DictReader((MODULE / "results/table_s6.csv").read_text().splitlines())
        actual = {row["dataset"]: tuple(f"{float(row[key]):.4f}" for key in diagnostic.TABLE_FIELDS[1:])
                  for row in rows}
        self.assertEqual(actual, published)


if __name__ == "__main__":
    unittest.main()
