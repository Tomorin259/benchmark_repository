#!/usr/bin/env python3
"""Verify registered external artifacts and create explicit per-package archives."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import tarfile

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    result = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def prepare_release(artifact_root, output_root, package=None, build=False, registry=None, mapping=None):
    with Path(registry or ROOT / 'artifacts.tsv').open(newline='') as handle:
        registered = {r['artifact_id']: r for r in csv.DictReader(handle, delimiter='\t') if r['status'] == 'external'}
    with Path(mapping or ROOT / 'data/release_map.tsv').open(newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    if {r['artifact_id'] for r in rows} != set(registered) or len(rows) != len(registered):
        raise ValueError('Release map must cover every external artifact exactly once')
    chosen = [r for r in rows if package is None or r['package_id'] == package]
    if not chosen:
        raise ValueError('Unknown or empty package')
    root, output = Path(artifact_root).resolve(), Path(output_root).resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError('Output directory must be empty')
    records, unresolved, verified_paths = [], [], {}
    for row in chosen:
        expected = registered[row['artifact_id']]['sha256']
        if not row['source_locator'] or not expected:
            unresolved.append({'artifact_id': row['artifact_id'], 'reason': row['source_status']})
            continue
        locator = row['source_locator']
        if locator.startswith('velocity:'):
            source_root = root
            relative = Path(locator.removeprefix('velocity:'))
        elif locator.startswith('repository:'):
            source_root = ROOT.resolve()
            relative = Path(locator.removeprefix('repository:'))
        else:
            raise ValueError(f'Unsupported locator: {locator}')
        path = (source_root / relative).resolve()
        if relative.is_absolute() or '..' in relative.parts or not path.is_relative_to(source_root) or not path.is_file():
            unresolved.append({'artifact_id': row['artifact_id'], 'reason': 'missing_source'})
            continue
        actual = digest(path)
        if actual != expected:
            unresolved.append({'artifact_id': row['artifact_id'], 'reason': 'hash_mismatch', 'expected': expected, 'actual': actual})
            continue
        member = Path(row['member_path'])
        if member.is_absolute() or '..' in member.parts or not row['member_path']:
            raise ValueError('Archive member must be a safe relative path')
        records.append({**row, 'sha256': expected, 'bytes': path.stat().st_size})
        verified_paths[row['artifact_id']] = path
    member_keys = [(r['package_id'], r['member_path']) for r in records]
    if len(set(member_keys)) != len(member_keys):
        raise ValueError('Duplicate archive member paths')
    if build and unresolved:
        raise ValueError(f'{len(unresolved)} unresolved sources; no archive may be published')
    output.mkdir(parents=True, exist_ok=True)
    archives = []
    if build:
        for identifier in sorted({r['package_id'] for r in records}):
            if not identifier or any(x in identifier for x in ('/', '\\', '..')):
                raise ValueError('Invalid package ID')
            target = output / (identifier + '.tar.gz')
            with tarfile.open(target, 'w:gz') as archive:
                for record in [r for r in records if r['package_id'] == identifier]:
                    archive.add(verified_paths[record['artifact_id']], arcname=record['member_path'], recursive=False)
            archives.append({'path': target.name, 'sha256': digest(target), 'bytes': target.stat().st_size})
    report = {'registered': len(chosen), 'verified': records, 'unresolved': unresolved,
              'archives': archives, 'published': False, 'download_uris': 'remain empty until actual publication'}
    (output / 'release_inventory.json').write_text(json.dumps(report, indent=2) + '\n')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--artifact-root', type=Path, required=True,
                        help='Root for velocity: sources; repository: sources use this checkout')
    parser.add_argument('--output-root', type=Path, required=True)
    parser.add_argument('--package')
    parser.add_argument('--build', action='store_true', help='Build only after every selected file passes hash verification')
    args = parser.parse_args()
    report = prepare_release(args.artifact_root, args.output_root, args.package, args.build)
    print(f"{len(report['verified'])}/{report['registered']} sources verified; {len(report['unresolved'])} unresolved; publication is not performed by this script")
