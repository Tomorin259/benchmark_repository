import csv
import hashlib
import json
import re
import subprocess
import sys
import tempfile
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
        with (ROOT / 'data/input_file_manifest.tsv').open() as handle:
            downloads = {r['artifact_id']: r['url'] for r in csv.DictReader(handle, delimiter='\t')}
        for row in self.read_rows():
            self.assertFalse(row["source_path"] if "source_path" in row else "")
            self.assertEqual(row["source_status"], "source_sha256_verified")
            if row["role"] == "reference_table":
                self.assertEqual(row["public_storage"], "upstream_only")
                self.assertEqual(row["public_url"], "-")
            else:
                self.assertEqual(row["public_storage"], "huggingface_file")
                self.assertEqual(row["public_url"], downloads[row['artifact_id']])
                self.assertTrue(row["public_url"].endswith(row["member_path"] + "?download=true"))
            self.assertFalse(re.search(r"/(?:mnt|Users|home)/|[A-Za-z]:\\", row["member_path"]))
            self.assertRegex(row["sha256"], r"^[0-9a-f]{64}$")
            self.assertGreater(int(row["bytes"]), 0)

    def test_individual_files_match_scientific_registry_without_original_tables(self):
        with (ROOT / 'data/input_file_manifest.tsv').open(newline='') as handle:
            files = list(csv.DictReader(handle, delimiter='\t'))
        with (ROOT / 'artifacts.tsv').open(newline='') as handle:
            registry = {r['artifact_id']: r for r in csv.DictReader(handle, delimiter='\t')}
        inputs = {r['artifact_id']: r for r in self.read_rows()}
        self.assertEqual(len(files), 25)
        self.assertEqual({r['artifact_id'] for r in files}, {
            r['artifact_id'] for r in inputs.values() if r['role'] != 'reference_table'})
        self.assertEqual({ext: sum(r['path'].endswith(ext) for r in files) for ext in ['.h5ad', '.loom', '.csv']},
                         {'.h5ad': 22, '.loom': 1, '.csv': 2})
        for entry in files:
            source = inputs[entry['artifact_id']]
            self.assertEqual(entry['sha256'], source['sha256'])
            self.assertEqual(entry['bytes'], source['bytes'])
            self.assertEqual(entry['path'], source['member_path'])
            self.assertEqual(entry['url'], source['public_url'])
            self.assertEqual(entry['url'], registry[entry['artifact_id']]['uri'])
            self.assertTrue(entry['url'].endswith(entry['path'] + '?download=true'))
            self.assertNotIn('aax3072_table-', entry['path'])

    def test_original_reference_tables_are_identified_but_not_mirrored(self):
        with (ROOT / 'data/upstream_references.tsv').open(newline='') as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        inputs = {r['artifact_id']: r for r in self.read_rows()}
        self.assertEqual({r['table'] for r in rows}, {'S1', 'S3', 'S4'})
        for row in rows:
            self.assertEqual(row['source_url'], 'https://doi.org/10.1126/science.aax3072')
            self.assertEqual(row['benchmark_csv_sha256'], inputs[row['artifact_id']]['sha256'])
            self.assertEqual(row['velocity_member_path'], inputs[row['artifact_id']]['member_path'])
            self.assertEqual(row['distribution'], 'upstream_only_not_mirrored')

    def test_reference_checker_rejects_missing_and_changed_bytes_without_writing(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = root / 'references.tsv'
            body = b'original reference fixture\n'
            with manifest.open('w', newline='') as handle:
                writer = csv.DictWriter(handle, delimiter='\t', fieldnames=[
                    'table', 'source_url', 'velocity_member_path', 'kinetics_member_path',
                    'benchmark_csv_bytes', 'benchmark_csv_sha256'])
                writer.writeheader()
                writer.writerow(dict(table='S1', source_url='https://doi.org/10.1126/science.aax3072',
                                     velocity_member_path='velocity/reference.csv',
                                     kinetics_member_path='kinetics/reference.csv',
                                     benchmark_csv_bytes=len(body),
                                     benchmark_csv_sha256=hashlib.sha256(body).hexdigest()))
            command = [sys.executable, str(ROOT / 'data/check_references.py'),
                       '--artifact-root', str(root), '--manifest', str(manifest)]
            missing = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(missing.returncode, 1)
            self.assertIn('MISSING S1', missing.stdout)
            self.assertFalse((root / 'velocity').exists())
            reference = root / 'kinetics/reference.csv'
            reference.parent.mkdir()
            reference.write_bytes(body)
            self.assertEqual(subprocess.run(command, capture_output=True).returncode, 0)
            reference.write_bytes(b'changed reference\n')
            changed = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(changed.returncode, 1)
            self.assertIn('MISMATCH S1', changed.stdout)
            self.assertEqual(reference.read_bytes(), b'changed reference\n')

    def test_expected_layer_contracts_are_visible(self):
        rows = {row["artifact_id"]: row for row in self.read_rows()}
        self.assertEqual(rows["input-labelled-sceu_rpe1"]["layers"], "sl;su;ul;uu")
        self.assertEqual(rows["input-labelled-sceu_organoid"]["layers"], "sl;su;ul;uu")
        self.assertEqual(set(rows["input-sim-only-murk-clean"]["layers"].split(';')),
                         {"alpha_true", "beta_true", "gamma_true", "spliced", "unspliced",
                          "clean_spliced", "clean_unspliced", "velocity_s_true", "velocity_u_true"})
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
