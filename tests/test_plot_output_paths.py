"""Plot exports stay in their requested directory and preserve source assets."""

import contextlib
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def load_script(relative):
    spec = importlib.util.spec_from_file_location(Path(relative).stem, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def hashes(paths):
    return {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


class PlotOutputTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cell = load_script("evaluation/cell/simulation/plot_simulation_cell_comparison.py")
        cls.cases = load_script("supplement/gene-kinetic/gene_kinetic_pilot.py")

    def assert_cli_error(self, module, arguments, message):
        stderr = io.StringIO()
        with patch.object(sys, "argv", [module.__file__, *arguments]), contextlib.redirect_stderr(stderr):
            with self.assertRaises(SystemExit) as error:
                module.main()
        self.assertEqual(error.exception.code, 2)
        self.assertIn(message, stderr.getvalue())

    def test_output_directory_is_required(self):
        for module in (self.cell, self.cases):
            with self.subTest(entry=module.__file__):
                self.assert_cli_error(module, [], "--output-dir")

    def test_protected_nonempty_and_file_destinations_are_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            occupied = Path(folder) / "occupied"
            occupied.mkdir()
            sentinel = occupied / "keep.txt"
            sentinel.write_bytes(b"previous export")
            regular_file = Path(folder) / "file"
            regular_file.write_bytes(b"input")
            for module in (self.cell, self.cases):
                for destination in (occupied, regular_file, ROOT / "evaluation/cell/simulation/figures"):
                    with self.subTest(entry=module.__file__, destination=destination):
                        self.assert_cli_error(module, ["--output-dir", str(destination)], "--output-dir")
            self.assertEqual(sentinel.read_bytes(), b"previous export")
            self.assertEqual(regular_file.read_bytes(), b"input")

    def test_cell_export_preserves_scores_and_published_bundle(self):
        module = self.cell
        protected = hashes([module.HERE / "results/cross_level_summary.csv",
                            *sorted((module.HERE / "figures").glob("simulation_*"))])
        expected = pd.read_csv(module.HERE / "results/cross_level_summary.csv")
        output_parent = ROOT / "outputs"
        output_parent.mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(dir=output_parent) as folder:
            output = Path(folder) / "panels"
            output.mkdir()
            arguments = ["--output-dir", str(output), "--cross-level-only"]
            with patch.object(sys, "argv", [module.__file__, *arguments]):
                module.main()
            self.assertEqual({p.name for p in output.iterdir()}, {
                "cross_level_summary.csv", "simulation_cross_level_cbdir.png",
                "simulation_cross_level_cbdir.pdf", "simulation_cross_level_cbdir.svg",
            })
            actual = pd.read_csv(output / "cross_level_summary.csv")
            pd.testing.assert_frame_equal(actual, expected, check_exact=False, rtol=1e-12, atol=1e-12)
            self.assertEqual(len(actual), 36)
            exported = hashes(output.iterdir())
            self.assert_cli_error(module, arguments, "must be an empty directory")
            self.assertEqual(hashes(exported), exported)
        self.assertEqual(hashes(protected), protected)

    def test_cases_resolve_relative_output_before_chdir_and_keep_diaph3_input(self):
        from evaluation.gene.config import REAL_CONFIGS

        module = self.cases
        names = ("dentate_diaph3", "dentate_ncald")
        selected_cases = {name: module.CASES[name] for name in names}
        diaph3_path = module.HERE / "results/diaph3_nipc_to_nbl1.csv"
        protected = hashes([diaph3_path, *sorted((module.HERE / "figures").iterdir())])
        with tempfile.TemporaryDirectory() as folder:
            workspace = Path(folder) / "invocation"
            workspace.mkdir()
            artifact = Path(folder) / "inputs/velocity"
            clusters = ["nIPC", "Nbl1", "Nbl2", "ImmGranule1", "ImmGranule2", "Granule"]
            values = np.arange(1, 73, dtype=float).reshape(36, 2) / 10
            data = ad.AnnData(
                X=values.copy(),
                obs=pd.DataFrame({"clusters": pd.Categorical(clusters * 6, categories=clusters)},
                                 index=[f"cell{i}" for i in range(36)]),
                var=pd.DataFrame(index=["Diaph3", "Ncald"]),
            )
            data.layers["Ms"] = values.copy()
            data.layers["Mu"] = values * 0.6 + 0.1
            for pair in module.VELOCITY_KEYS.values():
                for key in pair:
                    data.layers[key] = np.sin(values) * 0.1
            data.uns["clusters_colors"] = ["#112233", "#334455", "#556677",
                                           "#778899", "#99aabb", "#bbccdd"]
            inputs = []
            for relative in REAL_CONFIGS["DentateGyrus_lamanno"]["paths"].values():
                path = artifact / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                data.write_h5ad(path)
                inputs.append(path)
            original_inputs = hashes(inputs)
            arguments = ["--artifact-root", str(artifact), "--case", "all",
                         "--output-dir", "outputs/cases"]
            previous_cwd = Path.cwd()
            try:
                os.chdir(workspace)
                with patch.dict(module.CASES, selected_cases, clear=True), \
                        patch.object(sys, "argv", [module.__file__, *arguments]):
                    module.main()
            finally:
                os.chdir(previous_cwd)
            output = workspace / "outputs/cases"
            self.assertEqual({p.name for p in output.iterdir()}, {
                f"{name}_{suffix}" for name in names
                for suffix in ("local_phase.png", "local_phase.pdf", "manifest.json")
            })
            self.assertFalse((artifact.parent / "outputs").exists())
            manifest = json.loads((output / "dentate_diaph3_manifest.json").read_text())
            self.assertIn("supplement/gene-kinetic/results/diaph3_nipc_to_nbl1.csv",
                          manifest["edge_score_sources"])
            scores = pd.read_csv(diaph3_path)
            expected = scores.loc[scores.method.eq("cellDancer") & scores.source.eq("nIPC")
                                  & scores.target.eq("Nbl1") & scores.status.eq("ok"), "edge_score"].item()
            self.assertEqual(manifest["edge_scores"]["cellDancer"][0],
                             "nIPC→Nbl1: %+.3f" % expected)
            self.assertEqual(hashes(inputs), original_inputs)
            self.assert_cli_error(module, ["--output-dir", str(output)], "must be an empty directory")
        self.assertEqual(hashes(protected), protected)


if __name__ == "__main__":
    unittest.main()
