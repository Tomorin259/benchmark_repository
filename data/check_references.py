#!/usr/bin/env python3
"""Check user-supplied original references without downloading or rewriting them."""
import argparse
import csv
import hashlib
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--artifact-root', type=Path, required=True)
parser.add_argument('--manifest', type=Path, default=Path(__file__).with_name('upstream_references.tsv'))
args = parser.parse_args()
failed = False
with args.manifest.open(newline='') as handle:
    rows = list(csv.DictReader(handle, delimiter='\t'))
for row in rows:
    candidates = [args.artifact_root / row[key] for key in ('velocity_member_path', 'kinetics_member_path')]
    present = [p for p in candidates if p.is_file()]
    if not present:
        print(f"MISSING {row['table']}: {row['source_url']} -> {candidates[0]} or {candidates[1]}")
        failed = True
    for path in present:
        h = hashlib.sha256()
        with path.open('rb') as source:
            for block in iter(lambda: source.read(1024 * 1024), b''):
                h.update(block)
        match = path.stat().st_size == int(row['benchmark_csv_bytes']) and h.hexdigest() == row['benchmark_csv_sha256']
        print(f"{'OK' if match else 'MISMATCH'} {row['table']}: {path}")
        failed |= not match
raise SystemExit(1 if failed else 0)
