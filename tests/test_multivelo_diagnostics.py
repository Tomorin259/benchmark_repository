"""Lightweight source-data contracts for the post-hoc MultiVelo materials."""
import csv
import hashlib
import importlib.util
import json
import math
import re
import subprocess
import sys
import tempfile
import unittest
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
MODULE = REPO / "supplement/multivelo"
spec = importlib.util.spec_from_file_location("multivelo_summary", MODULE / "summarize_controls.py")
summary = importlib.util.module_from_spec(spec)
spec.loader.exec_module(summary)


class MultiVeloDiagnosticsTests(unittest.TestCase):
    def test_frozen_results_and_portable_provenance(self):
        manifest = json.loads((MODULE / "provenance.json").read_text())
        for name, expected in manifest["outputs"].items():
            self.assertEqual(hashlib.sha256((MODULE / name).read_bytes()).hexdigest(), expected, name)
        text = (MODULE / "provenance.json").read_text()
        self.assertNotIn("/mnt/", text)
        self.assertNotIn("/Users/", text)
        rows = summary.read(MODULE / "results/multivelo_refits.csv")
        self.assertEqual(len(rows), 250)
        self.assertEqual(Counter(r["condition"] for r in rows), {"refit_default": 70, "refit_complete": 70, "refit_no_output_smoothing": 70, "on_allow_repression": 40})
        self.assertEqual(len({(r["dataset"], r["gene"], r["condition"]) for r in rows}), 250)
        self.assertEqual(len(summary.read(MODULE / "results/scvelo_controls.csv")), 490)

    def test_default_reproduction_and_single_factor_checks(self):
        rows = summary.read(MODULE / "results/candidate_gate_baseline_reproduction.csv")
        self.assertEqual(len(rows), 40)
        self.assertTrue(all(float(r[k]) == 0 for r in rows for k in ["max_abs_velocity_error", "cbdir_error", "cossim_error"]))
        checks = json.loads((MODULE / "results/patch_validation.json").read_text())
        self.assertTrue(all(checks.values()))
        effects = summary.read(MODULE / "results/candidate_gate_effects.csv")
        self.assertEqual(len(effects), 2)
        self.assertTrue(all(int(r[k]) == 20 for r in effects for k in ["n_gene_pairs", "n_lower_cbdir", "n_lower_cossim", "n_lower_fit_loss"]))

    def test_summary_regeneration(self):
        with tempfile.TemporaryDirectory() as directory:
            subprocess.run([sys.executable, str(MODULE / "summarize_controls.py"), "--output-root", directory], check=True, capture_output=True)
            actual = summary.read(Path(directory) / "condition_summary.csv")
            expected = summary.read(MODULE / "results/condition_summary.csv")
            self.assertEqual(len(actual), 22)
            self.assertEqual(len(actual), len(expected))
            for rebuilt, frozen in zip(actual, expected):
                self.assertEqual(rebuilt.keys(), frozen.keys())
                for key in rebuilt:
                    if key.startswith(("mean_", "sd_")):
                        # statistics.stdev differs in the final binary digits
                        # across supported Python versions; source hashes stay exact.
                        self.assertTrue(math.isclose(float(rebuilt[key]), float(frozen[key]), rel_tol=0, abs_tol=1e-12), key)
                        self.assertEqual(f"{float(rebuilt[key]):.4f}", f"{float(frozen[key]):.4f}")
                    else:
                        self.assertEqual(rebuilt[key], frozen[key])
            by_key = {(r["dataset"], r["gene_class"], r["method"], r["condition"]): r for r in actual}
            self.assertEqual(f"{float(by_key['branching','branching','MultiVelo (RNA_only)','Gate open']['mean_cossim']):.4f}", "0.2477")
            self.assertEqual(f"{float(by_key['all','branching','MultiVelo (RNA_only)','Gate open']['mean_cossim']):.4f}", "0.2349")

    def test_phase_input_identity_and_no_pickle(self):
        try:
            import numpy as np
        except ImportError:
            self.skipTest("NumPy is required for array checks")
        with np.load(MODULE / "inputs/phase_case.npz", allow_pickle=False) as data:
            self.assertEqual(len(data["cell_ids"]), 2000)
            self.assertEqual(len(set(data["cell_ids"])), 2000)
            self.assertEqual(data["Ms"].shape, (2000, 1))
            self.assertEqual(data["Mu"].shape, (2000, 1))
            self.assertTrue(all(data[k].dtype.kind != "O" for k in data.files))
        scv = summary.read(MODULE / "results/scvelo_phase_assignment_case.csv")
        branch_b = [r for r in scv if r["bin"].startswith("BranchB_")]
        self.assertEqual(len(branch_b), 1000)
        self.assertEqual(sum(r["final_state"] == "repression" for r in branch_b), 540)

    def test_accepted_materials_and_plotted_scores(self):
        manifest = json.loads((MODULE / "materials_manifest.json").read_text())
        self.assertFalse(manifest["panel_letters"])
        self.assertTrue(manifest["input_hashes_unchanged"])
        self.assertFalse(manifest["primary_results_modified"])
        for name, expected in manifest["outputs"].items():
            self.assertEqual(hashlib.sha256((MODULE / name).read_bytes()).hexdigest(), expected, name)
        for name, expected in manifest["inputs"].items():
            self.assertEqual(hashlib.sha256((MODULE / name).read_bytes()).hexdigest(), expected, name)
        diagnostics = manifest["phase_diagnostics"]
        self.assertEqual(len(diagnostics), 3)
        self.assertTrue(all(row["sampled_cell_ids"] == diagnostics[0]["sampled_cell_ids"] for row in diagnostics))
        rows = summary.read(MODULE / "figures/scores_used.csv")
        self.assertEqual(len(rows), 160)
        sources = summary.read(MODULE / "results/multivelo_refits.csv") + summary.read(MODULE / "results/scvelo_controls.csv")
        by_key = {(r["dataset"], r["gene"], r.get("condition", r.get("variant"))): r for r in sources}
        for row in rows:
            source = by_key[row["dataset"], row["gene"], row["condition"]]
            self.assertAlmostEqual(float(row["cbdir"]), float(source["cbdir"]), places=14)
            self.assertAlmostEqual(float(row["cossim"]), float(source["cossim"]), places=14)

    def test_whole_gene_and_branch_score_scopes(self):
        cells = summary.read(MODULE / "results/candidate_gate_pilot_cell_states.csv")
        genes = summary.read(MODULE / "results/multivelo_refits.csv")
        branches = summary.read(MODULE / "results/candidate_gate_branch_scores.csv")
        for condition in ["refit_default", "on_allow_repression"]:
            whole = [r for r in cells if r["condition"] == condition]
            branch = [r for r in whole if r["branch"] == "BranchB"]
            self.assertEqual(len(whole), 2000)
            self.assertEqual(len(branch), 1000)
            for sample, source in [(whole, genes), (branch, branches)]:
                expected = next(r for r in source if r["dataset"] == "branching" and r["gene"] == "BoostALL_5" and r["condition"] == condition and (source is genes or r["branch"] == "BranchB"))
                actual = math.fsum(float(r["cossim"]) for r in sample) / len(sample)
                self.assertAlmostEqual(actual, float(expected["cossim"]), places=12)
            if condition == "on_allow_repression":
                self.assertEqual(sum(int(r["final_velocity_state"]) >= 2 for r in branch), 722)

    def test_manuscript_S15_and_S9_mapping(self):
        manifest = json.loads((REPO / "figures/suppfig15_composite_manifest.json").read_text())
        for item in [manifest["output"], manifest["assembly"], manifest["generator"], manifest["materials_manifest"], manifest["table"]] + manifest["panels"]:
            self.assertEqual(hashlib.sha256((REPO / item["path"]).read_bytes()).hexdigest(), item["sha256"])
        self.assertEqual({r["panel"] for r in manifest["panels"]}, {"A", "B"})
        self.assertEqual(manifest["table"]["latex_label"], "tab:suppMultiveloPhaseControls")
        table = (MODULE / "tables/tableS9.tex").read_text()
        rows = summary.read(MODULE / "results/condition_summary.csv")
        expected = []
        for method, condition, _ in summary.CONDITIONS[:4]:
            for dataset in ["branching", "all"]:
                row = next(r for r in rows if r["dataset"] == dataset and r["gene_class"] == "branching" and r["method"] == method and r["condition"] == condition)
                expected.extend(f"{float(row[key]):.4f}" for key in ["mean_cbdir", "mean_cossim"])
        self.assertEqual(re.findall(r"-?\d+\.\d{4}", table), expected)
        self.assertEqual(len(expected), 16)
        self.assertNotIn("20/20", table)
        self.assertNotIn("MURK", table)
        self.assertNotIn("smoothing", table)
        self.assertNotIn("Complete trajectory", table)
        self.assertEqual(table, summary.manuscript_table())
        material = json.loads((MODULE / "materials_manifest.json").read_text())
        self.assertIn("figures/paired_controls.pdf", material["outputs"])
        self.assertFalse(any("paired_branching" in path or "paired_all" in path for path in material["outputs"]))


if __name__ == "__main__":
    unittest.main()
