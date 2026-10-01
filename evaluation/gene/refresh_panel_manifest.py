#!/usr/bin/env python3
"""Check or refresh existing gene-panel hashes without rendering figures."""

import argparse
import csv
import hashlib
import io
from pathlib import Path


DEFAULT_MANIFEST = Path(__file__).resolve().parent / "figures/panels/assembly_manifest.tsv"
FIELDS = ["panel", "dataset", "asset", "path", "score_source", "sha256", "source_kind"]


def refresh_manifest(manifest: Path, write: bool = False) -> int:
    original = manifest.read_bytes()
    reader = csv.DictReader(io.StringIO(original.decode("utf-8")), delimiter="\t")
    if reader.fieldnames != FIELDS:
        raise ValueError("Unexpected gene-panel manifest columns")
    rows = list(reader)
    if not rows or any(set(row) != set(FIELDS) or None in row.values() for row in rows):
        raise ValueError("Empty or malformed gene-panel manifest")
    root = manifest.parent.resolve()
    hashes = {}
    changed = []
    identities = set()
    for row in rows:
        identity = tuple(row[field] for field in ("panel", "dataset", "asset", "path"))
        if identity in identities:
            raise ValueError(f"Duplicate panel reference: {identity}")
        identities.add(identity)
        relative = Path(row["path"])
        path = (root / relative).resolve()
        if not row["path"] or relative.is_absolute() or not path.is_relative_to(root):
            raise ValueError(f"Panel path must stay inside the manifest directory: {relative}")
        if not path.is_file():
            raise ValueError(f"Missing panel: {relative}")
        if row["path"] not in hashes:
            hashes[row["path"]] = hashlib.sha256(path.read_bytes()).hexdigest()
        if row["sha256"] != hashes[row["path"]]:
            changed.append(row["path"])
            row["sha256"] = hashes[row["path"]]
    if changed and not write:
        raise ValueError(f"Stale panel hashes in {len(changed)} references: {', '.join(changed)}")
    if changed:
        output = io.StringIO(newline="")
        newline = "\r\n" if b"\r\n" in original else "\n"
        writer = csv.DictWriter(output, fieldnames=FIELDS, delimiter="\t", lineterminator=newline)
        writer.writeheader()
        writer.writerows(rows)
        manifest.write_bytes(output.getvalue().encode("utf-8"))
    return len(changed)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--write", action="store_true", help="Refresh hashes; keep every other field")
    args = parser.parse_args()
    changed = refresh_manifest(args.manifest, args.write)
    print(f"Gene-panel manifest valid; {changed} hash references refreshed")
