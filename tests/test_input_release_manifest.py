import csv
import hashlib
import json
import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "data/input_release_manifest.tsv"
GENE_DIR = ROOT / "data/kinetics_reference_gene_usage"


class InputReleaseManifestTests(unittest.TestCase):
    def read_rows(self):
        with MANIFEST.open(newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def test_inventory_has_the_28_confirmed_source_files(self):
        rows = self.read_rows()
        self.assertEqual(len(rows), 28)
        self.assertEqual(len({row["artifact_id"] for row in rows}), 28)
        self.assertEqual(
            {row["role"] for row in rows},
            {"benchmark_main", "counts", "labelled_source", "fixed_simulation",
             "reference_table", "reference_gene_usage"},
        )
        self.assertEqual(sum(row["role"] == "benchmark_main" for row in rows), 14)
        self.assertEqual(sum(row["role"] == "counts" for row in rows), 4)
        self.assertEqual(sum(row["role"] == "labelled_source" for row in rows), 2)
        self.assertEqual(sum(row["role"] == "fixed_simulation" for row in rows), 3)
        self.assertEqual(sum(row["role"] == "reference_table" for row in rows), 3)
        self.assertEqual(sum(row["role"] == "reference_gene_usage" for row in rows), 2)

    def test_public_inventory_has_no_private_source_path_or_unverified_url(self):
        for row in self.read_rows():
            self.assertFalse(row["source_path"] if "source_path" in row else "")
            self.assertEqual(row["source_status"], "server_stat_sha256_verified_pending_publication")
            self.assertEqual(row["public_storage"], "zenodo_core_pending")
            self.assertEqual(row["public_url"], "-")
            self.assertFalse(re.search(r"/(?:mnt|Users|home)/|[A-Za-z]:\\", row["member_path"]))
            self.assertRegex(row["sha256"], r"^[0-9a-f]{64}$")
            self.assertGreater(int(row["bytes"]), 0)

    def test_expected_layer_contracts_are_visible(self):
        rows = {row["artifact_id"]: row for row in self.read_rows()}
        self.assertEqual(rows["input-labelled-sceu_rpe1"]["layers"], "sl;su;ul;uu")
        self.assertEqual(rows["input-labelled-sceu_organoid"]["layers"], "sl;su;ul;uu")
        self.assertEqual(rows["dataset-sim-only-murk"]["layers"],
                         "alpha_true;beta_true;gamma_true;spliced;unspliced;velocity_s_true;velocity_u_true")
        self.assertEqual(rows["dataset-mousebrain"]["n_obs"], "3365")
        self.assertEqual(rows["input-counts-mousebrain"]["n_obs"], "3365")

    def test_gene_usage_tables_match_their_manifest(self):
        metadata = json.loads((GENE_DIR / "manifest.json").read_text())
        for name, info in metadata["outputs"].items():
            path = GENE_DIR / name
            self.assertTrue(path.is_file())
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), info["sha256"])
            with path.open() as handle:
                self.assertEqual(sum(1 for _ in handle) - 1, info["rows"])


if __name__ == "__main__":
    unittest.main()
