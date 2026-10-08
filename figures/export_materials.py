#!/usr/bin/env python3
"""Verify canonical paper sources or copy a fresh PowerPoint/WPS asset bundle."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / 'figures/assembly/assembly_manifest.tsv'


def digest(path):
    result = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def export_materials(manifest=DEFAULT_MANIFEST, source_root=ROOT, output=None, figures=None):
    source_root = Path(source_root).resolve()
    with Path(manifest).open(newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    selected = [r for r in rows if not figures or r['figure'] in figures]
    if not selected:
        raise ValueError('No requested figures are registered')
    records, gaps, destinations = [], [], set()
    for row in selected:
        if row['status'] == 'missing_source':
            gaps.append({'figure': row['figure'], 'panel': row['panel'], 'note': row['note']})
            continue
        relative = Path(row['source_path'])
        source = (source_root / relative).resolve()
        if relative.is_absolute() or not source.is_relative_to(source_root) or not source.is_file():
            raise ValueError(f'Invalid or missing source: {relative}')
        if digest(source) != row['sha256']:
            raise ValueError(f'Source hash mismatch: {relative}')
        destination = Path(row['export_path'])
        if not row['export_path'] or destination.is_absolute() or '..' in destination.parts:
            raise ValueError(f'Invalid export path: {destination}')
        if destination.as_posix() in destinations:
            raise ValueError(f'Duplicate export destination: {destination}')
        destinations.add(destination.as_posix())
        records.append({'figure': row['figure'], 'panel': row['panel'],
                        'source': relative.as_posix(), 'export_path': destination.as_posix(),
                        'sha256': row['sha256']})
    if output is not None:
        output = Path(output).resolve()
        if output.is_relative_to(source_root) and not output.is_relative_to(source_root / 'outputs'):
            raise ValueError('Export outside the tracked repository or under outputs/')
        if output.exists() and any(output.iterdir()):
            raise ValueError('Export directory must be empty')
        output.mkdir(parents=True, exist_ok=True)
        for record in records:
            target = output / record['export_path']
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source_root / record['source'], target)
            if digest(target) != record['sha256']:
                raise ValueError(f'Export changed bytes: {target}')
        (output / 'bundle_manifest.json').write_text(json.dumps({
            'manifest_sha256': digest(Path(manifest)), 'files': records, 'source_gaps': gaps,
        }, indent=2) + '\n')
    return records, gaps


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument('--source-root', type=Path, default=ROOT)
    parser.add_argument('--figures', nargs='+')
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--verify-only', action='store_true')
    mode.add_argument('--output-root', type=Path)
    args = parser.parse_args()
    records, gaps = export_materials(args.manifest, args.source_root, args.output_root, args.figures)
    print(f'{len(records)} source files verified; {len(gaps)} explicit source gaps' +
          (f'; bundle exported to {args.output_root}' if args.output_root else ''))
