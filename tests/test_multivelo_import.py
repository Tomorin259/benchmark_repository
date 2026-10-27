"""Reject obsolete or incomplete bundles before changing formal results."""
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("multivelo_import", ROOT / "tools/archive_multivelo_results.py")
IMPORT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(IMPORT)


class MultiVeloImportTests(unittest.TestCase):
    def test_old_reference_cannot_replace_results(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            source, output = base / "source", base / "formal"
            source.mkdir(); output.mkdir()
            (output / "keep.csv").write_text("protected\n")
            (source / "provenance.json").write_text(json.dumps({"velocity_reference": {"sha256": "0" * 64}}))
            with self.assertRaisesRegex(ValueError, "current ODE reference"):
                IMPORT.import_results(source, output)
            self.assertEqual([p.name for p in output.iterdir()], ["keep.csv"])
            self.assertEqual((output / "keep.csv").read_text(), "protected\n")

    def test_subset_is_rejected_before_copying(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            source, output = base / "source", base / "formal"
            (source / "results").mkdir(parents=True); output.mkdir()
            truth = ROOT / "evaluation/gene/simulation/truth/manifest.json"
            (source / "provenance.json").write_text(json.dumps({"velocity_reference": {"sha256": IMPORT.sha(truth)}}))
            (source / "results/scvelo_controls.csv").write_text("dataset,gene,variant\nbranching,g,saved\n")
            with self.assertRaisesRegex(ValueError, "Incomplete"):
                IMPORT.import_results(source, output)
            self.assertEqual(list(output.iterdir()), [])


if __name__ == "__main__":
    unittest.main()
