#!/usr/bin/env python3
"""Check or refresh the seven-dataset finite common-gene CBDir comparison."""

import argparse
import csv
import math
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from evaluation.gene.config import METHOD_ORDER, REAL_CONFIGS

EXPECTED_COMMON_COUNTS = dict(zip(REAL_CONFIGS, (6, 6, 6, 7, 2, 4, 4)))


def common_gene_tables(source_path):
    with source_path.open(newline="") as handle:
        source = list(csv.DictReader(handle))
    finite, seen = {}, set()
    for row in source:
        if row["dataset"] not in REAL_CONFIGS:
            continue
        key = (row["dataset"], row["method"], row["gene"])
        if key in seen:
            raise ValueError(f"Duplicate real gene score: {key}")
        seen.add(key)
        value = float(row["new_cbdir"] or "nan")
        if row["status"] == "ok" and math.isfinite(value):
            finite[key] = value
    membership, dataset_rows = [], []
    for dataset in REAL_CONFIGS:
        by_method = [{g for d, m, g in finite if (d, m) == (dataset, method)}
                     for method in METHOD_ORDER]
        common = set.intersection(*by_method)
        if len(common) != EXPECTED_COMMON_COUNTS[dataset]:
            raise ValueError(f"Finite common-gene coverage changed: {dataset}, {len(common)}")
        membership.extend(dict(dataset=dataset, gene=gene) for gene in sorted(common))
        rows = []
        for method in METHOD_ORDER:
            values = [finite[dataset, method, gene] for gene in sorted(common)]
            rows.append(dict(dataset=dataset, method=method, n_common_genes=len(values),
                             mean=statistics.mean(values), sd=statistics.stdev(values)))
        for row in rows:
            row["rank"] = 1 + sum(other["mean"] > row["mean"] for other in rows)
        dataset_rows.extend(rows)
    macro = []
    for method in METHOD_ORDER:
        values = [row["mean"] for row in dataset_rows if row["method"] == method]
        macro.append(dict(method=method, n_datasets=len(values), mean=statistics.mean(values)))
    for row in macro:
        row["rank"] = 1 + sum(other["mean"] > row["mean"] for other in macro)
    return {"membership.csv": membership, "dataset_summary.csv": dataset_rows,
            "macro_summary.csv": macro}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=ROOT / "evaluation/gene/results/gene_cbdir_reference30_all.csv")
    parser.add_argument("--output-root", type=Path, default=ROOT / "evaluation/gene/results/common_genes")
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    for name, expected in common_gene_tables(args.source).items():
        path = args.output_root / name
        fields = list(expected[0])
        if args.write:
            path.parent.mkdir(parents=True, exist_ok=True)
            with path.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
                writer.writeheader()
                writer.writerows(expected)
        else:
            with path.open(newline="") as handle:
                reader = csv.DictReader(handle)
                actual = list(reader)
                if reader.fieldnames != fields or len(actual) != len(expected):
                    raise ValueError(f"Common-gene schema or coverage differs: {path}")
            for stored, computed in zip(actual, expected):
                for field in fields:
                    if field in ("mean", "sd"):
                        matches = math.isclose(float(stored[field]), computed[field], rel_tol=0, abs_tol=1e-12)
                    else:
                        matches = stored[field] == str(computed[field])
                    if not matches:
                        raise ValueError(f"Common-gene result differs: {path}, {computed}, {field}")
    print("Finite common genes: 6/6/6/7/2/4/4; 63 dataset means and 9 equal-dataset means")


if __name__ == "__main__":
    main()
