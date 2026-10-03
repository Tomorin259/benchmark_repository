import csv
import hashlib
import json
import tempfile
import unittest
from pathlib import Path


PLOTTING_ROOT = (
    Path(__file__).resolve().parents[1]
    / "evaluation"
    / "kinetics_recovery"
    / "plotting"
)


class KineticsPlottingContractTests(unittest.TestCase):
    def test_venn_two_line_titles_preserve_counts_and_membership(self):
        import pandas as pd
        from evaluation.kinetics_recovery.plotting.figure6 import plot_venn, venn2

        if venn2 is None:
            self.skipTest("matplotlib-venn is required for the Venn rendering check")

        root = PLOTTING_ROOT.parent
        source = root / "results/central_results"
        canonical = root / "figures/fig6/plot_data"
        config = json.loads((PLOTTING_ROOT / "figure6_config.json").read_text())
        with tempfile.TemporaryDirectory() as directory:
            for key, method, stem in (
                ("deepvelo_venn", "DeepVelo", "deepvelo"),
                ("cellDancer_venn", "cellDancer", "celldancer"),
            ):
                paths = config["inputs"]["rpe1"][key]
                output = Path(directory) / f"{stem}.svg"
                membership = Path(directory) / f"{stem}.csv"
                counts = plot_venn(
                    *(source / paths[name] for name in (
                        "reference_alpha", "reference_gamma",
                        "prediction_alpha", "prediction_gamma",
                    )),
                    output, method, membership,
                )
                pd.testing.assert_frame_equal(counts, pd.read_csv(
                    canonical / f"rpe1_{stem}_kinetic_classes.csv"))
                pd.testing.assert_frame_equal(pd.read_csv(membership), pd.read_csv(
                    canonical / f"rpe1_{stem}_kinetic_class_membership.csv"))
                svg = output.read_text()
                self.assertIn("Cooperative Genes", svg)
                self.assertIn("cosine", svg)
                self.assertNotIn("Macro-F1", svg)
                self.assertNotIn("cosine_similarity", svg)
                self.assertNotIn("Balanced accuracy", svg)

    def test_heatmap_means_preserve_scores_and_show_coverage(self):
        import pandas as pd
        from evaluation.kinetics_recovery.plotting.figure6 import plot_heatmap

        root = PLOTTING_ROOT.parent
        config = json.loads((PLOTTING_ROOT / "figure6_config.json").read_text())
        paths = config["inputs"]["rpe1"]["heatmap_metrics"]
        source = root / "results/central_results/gene_metrics.csv"
        source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "heatmap.svg"
            frames = plot_heatmap(
                {method: source for method in paths["Alpha"]},
                {method: source for method in paths["Gamma"]},
                config["heatmap_genes"], output,
            )
            for parameter, frame in zip(("alpha", "gamma"), frames):
                expected = pd.read_csv(
                    root / f"figures/fig6/plot_data/rpe1_{parameter}_heatmap.csv",
                    index_col=0,
                )
                pd.testing.assert_frame_equal(frame, expected, check_names=False)
            svg = output.read_text()
            for label in ("Dynamo (scEU)", "Mean (n)", "0.91", "(9)",
                          "0.42", "(6)"):
                self.assertIn(label, svg)
            for label in ("Rank / Method", "1. Dynamo", "Rows ranked by"):
                self.assertNotIn(label, svg)
        self.assertEqual(hashlib.sha256(source.read_bytes()).hexdigest(), source_hash)

    def test_trajectory_palette_separates_green_methods(self):
        from evaluation.cell.plots import BENCHMARK_COLORS
        from evaluation.kinetics_recovery.plotting.figure6 import TRAJECTORY_COLORS

        self.assertEqual(TRAJECTORY_COLORS["Dynamo(label data)"], "#126B36")
        self.assertEqual(TRAJECTORY_COLORS["cellDancer"], "#80C46F")
        for method in ("DeepVelo", "DeepKINET", "TSvelo"):
            self.assertEqual(TRAJECTORY_COLORS[method], BENCHMARK_COLORS[method])

    def test_organoid_boxplots_follow_benchmark_method_colors(self):
        from evaluation.cell.plots import BENCHMARK_COLORS
        from evaluation.kinetics_recovery.plotting.figure6 import BOXPLOT_METHOD_COLORS

        for method in ("DeepVelo", "cellDancer", "DeepKINET"):
            self.assertEqual(BOXPLOT_METHOD_COLORS[method], BENCHMARK_COLORS[method])
        self.assertEqual(
            BOXPLOT_METHOD_COLORS["Dynamo(label data)"],
            BENCHMARK_COLORS["Dynamo(RNA_seq)"],
        )

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

    def test_manual_stream_provenance_has_portable_input_locators(self):
        root = PLOTTING_ROOT.parents[2]
        manifests = {
            "evaluation/kinetics_recovery/figures/fig6/figure6E_organoid/"
            "dynamo_organoid_stream_title_removed.manifest.json":
                "promoted_canonical_figure6E_label_data_stream",
            "evaluation/kinetics_recovery/archive/fig5/panels/figure5E_organoid/"
            "dynamo_organoid_stream_title_removed.manifest.json":
                "review_only_not_manuscript_figure",
        }
        for relative, expected_status in manifests.items():
            manifest = json.loads((root / relative).read_text())
            self.assertEqual(manifest["status"], expected_status)
            for key in ("reference_h5ad", "dynamo_source_h5ad", "projection_csv"):
                locator = manifest[key]
                self.assertTrue(locator.startswith("artifact://velocity/"), locator)
                path = Path(locator.removeprefix("artifact://velocity/"))
                self.assertFalse(path.is_absolute(), locator)
                self.assertNotIn("..", path.parts, locator)

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
        self.assertEqual(
            set(manifest["stream_manifests"]),
            {"figure6A_rpe1", "figure6E_organoid"},
        )
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
        self.assertEqual(composite["export"], "validated_user_export_pdf")
        self.assertEqual(len(composite["composite_sha256"]), 64)
        self.assertEqual(len(composite["panels"]), 16)
        self.assertEqual(len({row["source_path"] for row in composite["panels"]}), 16)
        pending = set(composite["pending_panel_replacements"])
        self.assertEqual(
            composite["composite_status"],
            "pending_panel_replacements" if pending else "current",
        )
        for row in composite["panels"]:
            current_hash = hashlib.sha256((root / row["source_path"]).read_bytes()).hexdigest()
            if row["source_path"] in pending:
                self.assertEqual(current_hash, row["pending_source_sha256"])
                self.assertNotEqual(current_hash, row["source_sha256"])
            else:
                self.assertEqual(current_hash, row["source_sha256"])
        self.assertEqual(
            pending,
            {row["source_path"] for row in composite["panels"] if "pending_source_sha256" in row},
        )
        self.assertEqual(hashlib.sha256((root / composite["composite_path"]).read_bytes()).hexdigest(), composite["composite_sha256"])


if __name__ == "__main__":
    unittest.main()
