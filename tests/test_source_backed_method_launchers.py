import csv
import hashlib
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


def load_launcher(name):
    path = REPOSITORY_ROOT / "methods" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class SourceBackedLauncherTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tfvelo = load_launcher("tfvelo")
        cls.tsvelo = load_launcher("tsvelo")

    def test_both_launchers_cover_unique_fourteen_dataset_contracts(self):
        self.assertEqual(len(self.tfvelo.DATASETS), 14)
        self.assertEqual(len(self.tsvelo.DATASETS), 14)
        self.assertEqual(len({row[2] for row in self.tfvelo.DATASETS.values()}), 14)
        self.assertEqual(len({row[2] for row in self.tsvelo.DATASETS.values()}), 14)

    def test_official_fallbacks_are_explicitly_staged(self):
        self.assertEqual(self.tfvelo.STAGE_PATHS["Pancreas"].split("/")[-1], "endocrinogenesis_day15.h5ad")
        self.assertEqual(self.tfvelo.STAGE_PATHS["Gastrulation"].split("/")[-1], "erythroid_lineage.h5ad")
        self.assertEqual(self.tfvelo.DATASETS["ForebrainGlut"][3], "pca")
        self.assertNotIn(
            self.tsvelo.DATASETS["mousebrain"][2],
            {"gastrulation_erythroid", "pancreas", "pons", "10x_mouse_brain", "dentategyrus1", "larry"},
        )

    def test_tsvelo_mousebrain_contract_is_recorded(self):
        config = self.tsvelo.DATASETS["mousebrain"]
        self.assertTrue(config[3])
        self.assertIn("10X_multiome_mouse_brain.loom", config[0])
        self.assertEqual(self.tsvelo.normalize_cell_id("RNA_ATAC:AACCTTAAGTAGCGGGx"), "AACCTTAAGTAGCGGG-1")

    def test_source_patches_match_registered_revisions_and_hashes(self):
        methods = REPOSITORY_ROOT / "methods"
        sources = json.loads((methods / "source_manifest.json").read_text())["methods"]
        with (methods / "configurations.tsv").open() as handle:
            registry = {row["config_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
        self.assertEqual(set(sources), {"C10", "C13", "C14"})
        for config_id, source in sources.items():
            with self.subTest(config_id=config_id):
                self.assertRegex(source["revision"], r"^[0-9a-f]{40}$")
                self.assertIn(source["revision"], registry[config_id]["version"])
                self.assertEqual(source["repository"], registry[config_id]["official_source"])
                patch = methods / source["patch"]
                self.assertEqual(hashlib.sha256(patch.read_bytes()).hexdigest(), source["patch_sha256"])
                self.assertTrue(source["files"])
                for fingerprint in source["files"].values():
                    self.assertRegex(fingerprint, r"^[0-9a-f]{64}$")

    def test_static_source_and_input_validation(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            tf_source = root / "tfvelo"
            for relative in (
                "TFvelo_run_demo.py",
                "TFvelo_analysis_demo.py",
                "TFvelo/__init__.py",
                "ENCODE/processed",
                "ChEA/ChEA_2016.txt",
            ):
                path = tf_source / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.touch()
            tf_data = root / "tf_data"
            for row in self.tfvelo.DATASETS.values():
                path = tf_data / row[0]
                path.parent.mkdir(parents=True, exist_ok=True)
                path.touch()
            self.tfvelo.check_paths(tf_source, tf_data, self.tfvelo.DATASETS_TO_RUN)

            ts_source = root / "tsvelo"
            for relative in (
                "TSvelo_run.py",
                "TSvelo/TSvelo_pp.py",
                "TSvelo/TSvelo_branch.py",
                "TSvelo/TSvelo_concate.py",
                "ENCODE/processed",
                "ChEA/ChEA_2016.txt",
            ):
                path = ts_source / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.touch()
            ts_data = root / "ts_data"
            for row in self.tsvelo.DATASETS.values():
                path = ts_data / row[0]
                path.parent.mkdir(parents=True, exist_ok=True)
                path.touch()
            self.tsvelo.check_paths(ts_source, ts_data, self.tsvelo.DATASETS_TO_RUN)

            for name, launcher, source, data in (
                ("tfvelo", self.tfvelo, tf_source, tf_data),
                ("tsvelo", self.tsvelo, ts_source, ts_data),
            ):
                with self.subTest(launcher=name):
                    for dataset, row in launcher.DATASETS.items():
                        if dataset != "Pancreas":
                            (data / row[0]).unlink()
                    launcher.check_paths(source, data, ["Pancreas"])
                    with self.assertRaises(FileNotFoundError):
                        launcher.check_paths(source, data, launcher.DATASETS_TO_RUN)
                    result = subprocess.run(
                        [sys.executable, "-B", str(REPOSITORY_ROOT / "methods" / f"{name}.py"),
                         "--dataset", "Pancreas", "--source-root", str(source),
                         "--data-root", str(data), "--check-only"],
                        capture_output=True, text=True,
                    )
                    self.assertEqual(result.returncode, 0, result.stderr)
                    self.assertIn("configuration valid: 1 datasets", result.stdout)
                    (data / launcher.DATASETS["Pancreas"][0]).unlink()
                    with self.assertRaises(FileNotFoundError):
                        launcher.check_paths(source, data, ["Pancreas"])


if __name__ == "__main__":
    unittest.main()
