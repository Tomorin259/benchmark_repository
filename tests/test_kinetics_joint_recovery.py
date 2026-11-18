import csv
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest

from evaluation.kinetics_recovery.joint_recovery import compare


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "evaluation/kinetics_recovery/results/joint_recovery"


def score_rows(dataset, branch, gene, values):
    return [
        {"Dataset": dataset, "Branch": branch, "Gene": gene, "Method": method,
         "Parameter": parameter, "Spearman_Corr": str(value)}
        for (method, parameter), value in zip(
            (("cellDancer", "alpha"), ("cellDancer", "gamma"),
             ("DeepVelo", "alpha"), ("DeepVelo", "gamma")),
            values,
        )
    ]


class JointRecoveryTests(unittest.TestCase):
    def test_intersection_and_equal_branch_weights(self):
        rows = score_rows("rpe1", "", "A", [0, 0, 0, 0])
        rows += score_rows("organoid", "1", "A", [0.2, 1, 0.6, 0])
        rows += score_rows("organoid", "1", "B", [0, 0, 0, "inf"])
        rows += score_rows("organoid", "1", "C", [0, "nan", 0, 0])
        rows += score_rows("organoid", "1", "D", [0, 0, 0])
        for gene in ("B", "C", "D"):
            rows += score_rows("organoid", "2", gene, [0, 0.2, 0.2, 0.5])
        scores, results = compare(rows)
        self.assertEqual(len(scores), 10)
        self.assertEqual({row[2] for row in scores if row[1] == "1"}, {"A"})
        self.assertEqual([c["common_genes"] for c in results["organoid"]["components"]], [1, 3])
        self.assertEqual(results["rpe1"]["front"], ["cellDancer", "DeepVelo"])
        self.assertEqual(results["organoid"]["front"], ["cellDancer", "DeepVelo"])
        self.assertAlmostEqual(results["organoid"]["means"]["cellDancer"]["gamma"], 0.6)
        self.assertAlmostEqual(results["organoid"]["means"]["DeepVelo"]["alpha"], 0.4)
        with self.assertRaisesRegex(ValueError, "Duplicate score"):
            compare(rows + rows[:1])

    def test_public_outputs_rebuild_from_canonical_scores(self):
        with TemporaryDirectory() as directory:
            output = Path(directory) / "comparison"
            command = [
                sys.executable, "-m", "evaluation.kinetics_recovery.joint_recovery",
                "--output", str(output),
            ]
            subprocess.run(command, cwd=ROOT, check=True)
            for name in ("common_gene_scores.csv", "common_gene_comparison.json"):
                self.assertEqual((output / name).read_bytes(), (RESULTS / name).read_bytes())
            with (output / "common_gene_scores.csv").open(newline="", encoding="utf-8") as handle:
                self.assertEqual(sum(1 for _ in csv.DictReader(handle)), 1578)
            results = json.loads((output / "common_gene_comparison.json").read_text())["results"]
            for system, expected in (("rpe1", [249]), ("organoid", [270, 270])):
                self.assertEqual([c["common_genes"] for c in results[system]["components"]], expected)
                self.assertEqual(results[system]["front"], ["DeepVelo"])
            rerun = subprocess.run(command, cwd=ROOT, capture_output=True)
            self.assertNotEqual(rerun.returncode, 0)
            self.assertIn(b"output directory must be new", rerun.stderr)


if __name__ == "__main__":
    unittest.main()
