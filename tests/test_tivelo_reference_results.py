"""Complete retained TIvelo controls and their shared truth identities."""
import csv
import hashlib
import json
import statistics
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE = ROOT / "supplement/tivelo_gene"


def read(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


class TIveloReferenceResultsTests(unittest.TestCase):
    def test_complete_root_pairs_and_summaries(self):
        for folder in ["root_stage_intervention", "root_stage_intervention_seed20260927"]:
            rows = read(MODULE / folder / "root_stage_gene_scores.csv")
            self.assertEqual(len(rows), 210)
            self.assertEqual(len({(r["scenario"], r["variant"], r["gene"]) for r in rows}), 210)
            summaries = read(MODULE / folder / "root_stage_summary.csv")
            self.assertEqual(len(summaries), 9)
            for summary in summaries:
                group = [r for r in rows if (r["scenario"], r["variant"]) == (summary["scenario"], summary["variant"])]
                self.assertEqual(len(group), int(summary["n_genes"]))
                for metric, field in [("cossim", "mean_cossim"), ("cbdir_reference30", "mean_cbdir_reference30")]:
                    self.assertAlmostEqual(statistics.mean(float(r[metric]) for r in group), float(summary[field]), places=12)

    def test_references_and_output_hashes(self):
        reference = ROOT / "evaluation/gene/simulation/truth/manifest.json"
        expected = hashlib.sha256(reference.read_bytes()).hexdigest()
        for folder, name in [("root_stage_intervention", "manifest.json"),
                             ("root_stage_intervention_seed20260927", "manifest.json"),
                             ("results", "root_intervention_manifest.json"),
                             ("results", "source_manifest.json"), ("stage_audit", "manifest.json"),
                             ("upstream_scvelo_audit", "manifest.json")]:
            manifest = json.loads((MODULE / folder / name).read_text())
            self.assertEqual(manifest["velocity_reference"]["sha256"], expected)
            for output, checksum in manifest["score_outputs_sha256"].items():
                self.assertEqual(hashlib.sha256((MODULE / folder / output).read_bytes()).hexdigest(), checksum)
        self.assertEqual(len(read(MODULE / "results/root_intervention_gene_scores.csv")), 720)
        self.assertEqual(len(read(MODULE / "stage_audit/stage_final_gene_scores.csv")), 70)

    def test_plotted_scores_and_main_cosine_replay_agree(self):
        roots = {(r["scenario"], r["variant"], r["gene"]): r
                 for r in read(MODULE / "root_stage_intervention/root_stage_gene_scores.csv")}
        plotted = read(MODULE / "figures/simulation_plotted_pairs.csv")
        self.assertEqual(len(plotted), 240)
        for row in plotted:
            source = roots[row["scenario"], row["variant"], row["gene"]]
            self.assertEqual(float(source[row["metric"]]), float(row["value"]))
        checks = read(MODULE / "results/cosine_recalculation_check.csv")
        self.assertEqual(len(checks), 70)
        self.assertLess(max(float(r["absolute_difference"]) for r in checks), 1e-12)


if __name__ == "__main__":
    unittest.main()
