"""Source identity, safe destinations and exact bytes in material exports."""
import csv
import hashlib
from pathlib import Path
import tempfile
import unittest

from figures.export_materials import export_materials


class MaterialExportTests(unittest.TestCase):
    def test_exact_export_and_tampered_source_rejection(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder) / 'source'
            root.mkdir()
            panel = root / 'panel.png'
            panel.write_bytes(b'canonical scientific panel')
            manifest = root / 'assembly.tsv'
            row = dict(figure='fig6', panel='C', role='heatmap', status='referenced_source',
                       source_path='panel.png', export_path='fig6/C.png',
                       sha256=hashlib.sha256(panel.read_bytes()).hexdigest(), note='')
            with manifest.open('w', newline='') as handle:
                writer = csv.DictWriter(handle, fieldnames=list(row), delimiter='\t')
                writer.writeheader(); writer.writerow(row)
            records, gaps = export_materials(manifest, root)
            self.assertEqual(len(records), 1); self.assertFalse(gaps)
            output = Path(folder) / 'bundle'
            export_materials(manifest, root, output)
            self.assertEqual((output / 'fig6/C.png').read_bytes(), panel.read_bytes())
            with self.assertRaisesRegex(ValueError, 'must be empty'):
                export_materials(manifest, root, output)
            panel.write_bytes(b'changed panel')
            with self.assertRaisesRegex(ValueError, 'hash mismatch'):
                export_materials(manifest, root, Path(folder) / 'new-bundle')
            self.assertFalse((Path(folder) / 'new-bundle').exists())


if __name__ == '__main__':
    unittest.main()
