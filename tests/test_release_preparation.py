import csv
import hashlib
import importlib.util
from pathlib import Path
import tarfile
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('prepare_release', ROOT / 'data/prepare_release.py')
release = importlib.util.module_from_spec(spec)
spec.loader.exec_module(release)


class ReleasePreparationTests(unittest.TestCase):
    def test_verified_bytes_are_packaged_and_changed_inputs_block_build(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            payload = b'frozen scientific input\n'
            (root / 'input.npz').write_bytes(payload)
            registry = root / 'registry.tsv'
            mapping = root / 'map.tsv'
            registry.write_text('artifact_id\tstatus\tsha256\ninput\texternal\t' + hashlib.sha256(payload).hexdigest() + '\n')
            mapping.write_text('artifact_id\tpackage_id\tsource_locator\tmember_path\tsource_status\ninput\ttoy\tvelocity:input.npz\tinputs/input.npz\trecovered\n')
            report = release.prepare_release(root, root / 'verified', build=True, registry=registry, mapping=mapping)
            self.assertEqual(len(report['verified']), 1)
            self.assertFalse(report['published'])
            with tarfile.open(root / 'verified/toy.tar.gz') as archive:
                self.assertEqual(archive.extractfile('inputs/input.npz').read(), payload)
            (root / 'input.npz').write_bytes(b'changed')
            with self.assertRaisesRegex(ValueError, 'unresolved sources'):
                release.prepare_release(root, root / 'blocked', build=True, registry=registry, mapping=mapping)
            self.assertFalse((root / 'blocked').exists())
            report = release.prepare_release(root, root / 'audit', registry=registry, mapping=mapping)
            self.assertEqual(report['unresolved'][0]['reason'], 'hash_mismatch')
            mapping.write_text(mapping.read_text().replace('inputs/input.npz', '../escape.npz'))
            (root / 'input.npz').write_bytes(payload)
            with self.assertRaisesRegex(ValueError, 'safe relative'):
                release.prepare_release(root, root / 'escape', build=True, registry=registry, mapping=mapping)

    def test_release_index_covers_external_registry_without_inventing_inputs(self):
        with (ROOT / 'artifacts.tsv').open(newline='') as handle:
            registry = {r['artifact_id']: r for r in csv.DictReader(handle, delimiter='\t') if r['status'] == 'external'}
        with (ROOT / 'data/release_map.tsv').open(newline='') as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        self.assertEqual(len(rows), len(registry))
        self.assertEqual({r['artifact_id'] for r in rows}, set(registry))
        for row in rows:
            if not registry[row['artifact_id']]['sha256']:
                self.assertEqual(row['source_locator'], '')
                self.assertEqual(row['source_status'], 'input_identity_not_frozen')
            else:
                self.assertTrue(row['source_locator'].startswith(('velocity:', 'repository:')))

    def test_repository_sources_are_packaged_from_the_checkout(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            checkout, artifacts = root / 'checkout', root / 'artifacts'
            checkout.mkdir()
            artifacts.mkdir()
            payload = b'registered repository output\n'
            (checkout / 'output.h5ad').write_bytes(payload)
            (artifacts / 'output.h5ad').write_bytes(b'different artifact-root file\n')
            registry, mapping = root / 'registry.tsv', root / 'map.tsv'
            registry.write_text('artifact_id\tstatus\tsha256\noutput\texternal\t' + hashlib.sha256(payload).hexdigest() + '\n')
            mapping.write_text('artifact_id\tpackage_id\tsource_locator\tmember_path\tsource_status\noutput\ttoy\trepository:output.h5ad\tmultivelo/output.h5ad\tverified\n')
            with patch.object(release, 'ROOT', checkout):
                report = release.prepare_release(artifacts, root / 'packed', build=True,
                                                 registry=registry, mapping=mapping)
            self.assertEqual(len(report['verified']), 1)
            with tarfile.open(root / 'packed/toy.tar.gz') as archive:
                self.assertEqual(archive.extractfile('multivelo/output.h5ad').read(), payload)


if __name__ == '__main__':
    unittest.main()
