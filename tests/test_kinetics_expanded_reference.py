import hashlib
from io import BytesIO
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from evaluation.kinetics_recovery.expanded_reference import (
    CANONICAL, RESULTS, candidate_genes, load_scores, summaries,
)


class ExpandedReferenceTests(unittest.TestCase):
    def test_scope_and_isolation(self):
        scope, scores = load_scores()
        self.assertEqual(len(scope), 829)
        self.assertFalse(scope.duplicated(["Dataset", "Gene"]).any())
        self.assertEqual(len(candidate_genes(scope, "rpe1")), 250)
        self.assertEqual(len(candidate_genes(scope, "organoid")), 270)
        org = scope[scope.Dataset.eq("organoid")]
        self.assertEqual(int(org.Source_selected.sum()), 295)
        self.assertEqual(int(org.Exclusion_reason.eq("housekeeping_control").sum()), 6)
        self.assertEqual(len(scores), 4748)
        self.assertEqual(int(np.isfinite(scores.Spearman_Corr).sum()), 4567)
        canonical = pd.read_csv(CANONICAL / "gene_metrics.csv")
        self.assertEqual(len(canonical), 4736)
        self.assertFalse((canonical.Method.eq("TSvelo") & canonical.Dataset.eq("organoid")).any())

    def test_tsvelo_mapping_and_undefined_scores(self):
        frame = pd.read_csv(RESULTS / "tsvelo_organoid_gene_metrics.csv")
        mapping = pd.read_csv(RESULTS / "tsvelo_organoid_gene_mapping.csv")
        self.assertEqual(len(mapping), 37)
        self.assertTrue(mapping.original_id_match.all())
        self.assertTrue(mapping.gene_id.str.startswith("ENSMUSG").all())
        self.assertEqual(set(mapping.reference_symbol), set(frame.Gene))
        self.assertTrue(frame.Parameter.eq("alpha").all())
        for branch, expected_mean, pairs in ((1, 0.09711947450027074, 196),
                                             (2, 0.2228012553949407, 187)):
            rows = frame[frame.Branch.eq(branch)]
            self.assertEqual(len(rows), 37)
            self.assertEqual(int(rows.Spearman_Corr.notna().sum()), 20)
            self.assertEqual(int(rows.Undefined_reason.eq("constant_prediction").sum()), 17)
            self.assertTrue(rows.Valid_Samples.eq(pairs).all())
            self.assertAlmostEqual(rows.Spearman_Corr.mean(), expected_mean, places=12)

    def test_summary_rebuild_and_common_membership(self):
        scope, scores = load_scores()
        rebuilt = summaries(scope, scores)
        for name, actual in rebuilt.items():
            expected = pd.read_csv(RESULTS / f"{name}.csv")
            pd.testing.assert_frame_equal(actual, expected, check_dtype=False, atol=1e-12, rtol=1e-12)
        coverage = rebuilt["coverage_summary"]
        self.assertEqual(len(coverage), 21)
        self.assertTrue((coverage.Valid_N + coverage.Undefined_N == coverage.Present_N).all())
        self.assertTrue((coverage.Present_N + coverage.Unavailable_N == coverage.Candidate_N).all())
        common = rebuilt["common_gene_summary"]
        self.assertEqual(len(common), 54)
        members = rebuilt["common_gene_membership"]
        for (panel, group), rows in common.groupby(["Panel", "Comparison_scope"]):
            chosen = members[members.Panel.eq(panel) & members.Comparison_scope.eq(group)]
            self.assertEqual(len(chosen), rows.Common_genes_N.iloc[0])
            self.assertTrue(rows.Method_valid_N.eq(len(chosen)).all())

    def test_manifest_hashes_and_portable_external_inputs(self):
        root = RESULTS.parents[3]
        manifest = json.loads((RESULTS / "manifest.json").read_text())
        for path, expected in manifest["outputs_sha256"].items():
            self.assertEqual(hashlib.sha256((root / path).read_bytes()).hexdigest(), expected, path)
        for record in manifest["external_inputs"].values():
            self.assertFalse(Path(record["path"]).is_absolute())
            self.assertNotIn("..", Path(record["path"]).parts)
        self.assertFalse(manifest["policies"]["canonical_scores_modified"])
        self.assertFalse(manifest["policies"]["figure7_aggregation_modified"])

    def test_trajectory_reserves_last_gene_and_complete_legend(self):
        import matplotlib.pyplot as plt
        from evaluation.kinetics_recovery.plotting.figure6 import plot_trajectory
        genes = [f"G{i}" for i in range(10)]
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            reference = directory / "reference.csv"
            method = directory / "method.csv"
            pd.DataFrame([[1, 2, 3]]*10, index=genes).to_csv(reference)
            # The first gene lacks this method: legend must not rely on it.
            pd.DataFrame([[3, 2, 1]]*9, index=genes[1:]).to_csv(method)
            with patch.object(plt, "close"):
                plotted = plot_trajectory("alpha", reference, {"TSvelo": method}, genes,
                                          BytesIO(), complete_legend=True, output_format="svg")
                figure = plt.gcf()
            self.assertEqual([axis.get_title() for axis in figure.axes if axis.axison], genes)
            self.assertEqual(set(plotted.Gene), set(genes))
            self.assertEqual([text.get_text() for text in figure.axes[-1].get_legend().texts],
                             ["Ground Truth", "TSvelo"])
            plt.close(figure)

    def test_boxplot_n_counts_only_finite_genes(self):
        import matplotlib.pyplot as plt
        from evaluation.kinetics_recovery.plotting.figure6 import plot_branch
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            source = directory / "scores.csv"
            pd.DataFrame([
                {"Method": "DeepVelo", "Dataset": "rpe1", "Branch": np.nan,
                 "Parameter": parameter, "Gene": gene, "Spearman_Corr": value}
                for parameter in ("alpha", "gamma") for gene, value in (("G1", .4), ("G2", np.nan))
            ]).to_csv(source, index=False)
            with patch.object(plt, "close"):
                plotted = plot_branch(None, {"DeepVelo": {"Alpha": source, "Gamma": source, "color": "red"}},
                                      ["G1", "G2"], directory / "plot.svg", dataset="rpe1",
                                      system_label="RPE1", show_n=True, dpi=30, bbox_inches=None)
                figure = plt.gcf()
            self.assertEqual(len(plotted), 4)
            self.assertEqual(list(figure.get_size_inches()), [18, 8])
            for axis in figure.axes:
                self.assertEqual(axis.get_xticklabels()[0].get_text(), "DeepVelo\n(n = 1)")
            plt.close(figure)


if __name__ == "__main__":
    unittest.main()
