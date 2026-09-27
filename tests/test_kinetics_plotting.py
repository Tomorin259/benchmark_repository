import csv
import hashlib
import json
import unittest
from pathlib import Path


PLOTTING_ROOT = (
    Path(__file__).resolve().parents[1]
    / "evaluation"
    / "kinetics_recovery"
    / "plotting"
)


class KineticsPlottingContractTests(unittest.TestCase):
    def test_stream_config_has_explicit_method_keys(self):
        config = json.loads((PLOTTING_ROOT / "stream_config.json").read_text())
        specs = config["organoid"]
        self.assertEqual([spec["method"] for spec in specs], ["cellDancer", "DeepKINET", "DeepVelo"])
        self.assertEqual(
            [spec["output_stem"] for spec in specs],
            ["celldancer", "DeepKINET", "deepvelo"],
        )
        self.assertEqual(
            [spec["velocity_layer"] for spec in specs],
            ["velocity", "DeepKINET_velocity", "velocity"],
        )
        self.assertEqual([spec["method"] for spec in config["rpe1_manual"]], ["Dynamo"])
        self.assertEqual([spec["method"] for spec in config["organoid_manual"]], ["Dynamo"])
        rpe1_specs = config["rpe1"]
        self.assertEqual([spec["method"] for spec in rpe1_specs],
                         ["cellDancer", "DeepVelo", "DeepKINET", "TSvelo"])
        self.assertEqual([spec["output_stem"] for spec in rpe1_specs],
                         ["celldancer", "deepvelo", "deepkinet", "tsvelo"])
        self.assertEqual([spec["velocity_layer"] for spec in rpe1_specs],
                         ["velocity", "velocity", "DeepKINET_velocity", "velocity"])
        for panel in ("rpe1", "organoid"):
            self.assertTrue(all(spec["result"].startswith("results/")
                                for spec in config[panel]))

    def test_stream_paths_are_relative(self):
        config = json.loads((PLOTTING_ROOT / "stream_config.json").read_text())
        for panel in config.values():
            for spec in panel:
                path = Path(spec.get("result") or spec.get("asset"))
                self.assertFalse(path.is_absolute())
                self.assertNotIn("..", path.parts)

    def test_kinetic_class_panels_use_mapped_trajectories(self):
        config = json.loads(
            (PLOTTING_ROOT / "figure6_config.json").read_text()
        )
        venn = config["inputs"]["rpe1"]["cellDancer_venn"]
        self.assertEqual(
            set(venn),
            {"reference_alpha", "reference_gamma", "prediction_alpha", "prediction_gamma"},
        )
        self.assertTrue(
            all("static" not in Path(path).parts for path in venn.values())
        )

    def test_promoted_figure_sources_and_assembly_match(self):
        root = PLOTTING_ROOT.parents[2]
        manifest = json.loads((PLOTTING_ROOT.parent / "results/central_results/manifest.json").read_text())
        for relative, expected in manifest["outputs_sha256"].items():
            self.assertNotIn("fig5_central_results", relative)
            self.assertEqual(hashlib.sha256((root / relative).read_bytes()).hexdigest(), expected, relative)
        with (root / "figures/assembly/assembly_manifest.tsv").open() as handle:
            rows = [row for row in csv.DictReader(handle, delimiter="\t")
                    if row["figure"] == "fig6" and row["status"] == "materialized"]
        self.assertEqual(len(rows), 16)
        expected_panel = {
            "organoid_branch1_boxplot.png": "F",
            "organoid_branch2_boxplot.png": "F",
            "rpe1_alpha_gamma_heatmap.png": "C",
            "rpe1_alpha_trajectories.png": "B",
            "rpe1_gamma_trajectories.png": "B",
            "rpe1_celldancer_kinetic_classes.png": "D",
            "rpe1_deepvelo_kinetic_classes.png": "D",
        }
        for row in rows:
            if Path(row["source_path"]).name in expected_panel:
                self.assertEqual(row["panel"], expected_panel[Path(row["source_path"]).name])
            source = (root / row["source_path"]).read_bytes()
            self.assertEqual(source, (root / "figures/assembly" / row["materialized_path"]).read_bytes())
            self.assertEqual(hashlib.sha256(source).hexdigest(), row["sha256"])
        composite = json.loads((root / "figures/fig6_composite_manifest.json").read_text())
        self.assertEqual(composite["paper_figure"], "fig6")
        self.assertEqual(composite["export"], "online_authoritative_pdf")
        self.assertEqual(len(composite["composite_sha256"]), 64)
        self.assertEqual(len(composite["panels"]), 16)
        self.assertEqual(len({row["source_path"] for row in composite["panels"]}), 16)
        self.assertEqual(composite["composite_status"], "current")
        self.assertFalse(composite["pending_panel_replacements"])
        for row in composite["panels"]:
            current_hash = hashlib.sha256((root / row["source_path"]).read_bytes()).hexdigest()
            self.assertEqual(current_hash, row["source_sha256"])
        self.assertEqual(hashlib.sha256((root / composite["composite_path"]).read_bytes()).hexdigest(), composite["composite_sha256"])


if __name__ == "__main__":
    unittest.main()
