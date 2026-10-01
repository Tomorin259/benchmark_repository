#!/usr/bin/env python3
"""Validate and export the simulated cell-level supplementary table.

The score CSV is the frozen 27-row source used by Supplemental Table S7.
The script does not rerun velocity methods or regenerate the user-exported
Supplemental Figure S12 PDF.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

METHODS = [
    "scVelo_dyn", "MultiVelo (RNA_only)", "UniTVelo", "TIvelo", "veloVI",
    "cellDancer", "DeepVelo", "DeepVelo_scv", "GraphVelo_scv",
]
DATASETS = ["Normal + MURK", "Normal + branching", "Hybrid"]
METRICS = ["Consistency Score", "CBDir", "ICCoh", "TransCosine", "BATC"]


def read_scores(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    expected = {"Method", "Dataset", *METRICS}
    if set(rows[0]) != expected:
        raise ValueError(f"Unexpected columns: {list(rows[0])}")
    if len(rows) != len(METHODS) * len(DATASETS):
        raise ValueError(f"Expected 27 rows, found {len(rows)}")
    keys = {(row["Method"], row["Dataset"]) for row in rows}
    if len(keys) != len(rows):
        raise ValueError("Duplicate method-dataset rows")
    if {row["Method"] for row in rows} != set(METHODS):
        raise ValueError("Method set does not match the frozen S12/S7 contract")
    if {row["Dataset"] for row in rows} != set(DATASETS):
        raise ValueError("Dataset set does not match the frozen S12/S7 contract")
    for row in rows:
        for metric in METRICS:
            value = row[metric]
            if len(value.split(".")[-1]) != 4:
                raise ValueError(f"Non-four-decimal value: {value}")
            float(value)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scores", type=Path, default=Path(__file__).parents[1] / "results/simulation/simulation_cell_scores_4dp.csv")
    parser.add_argument("--latex", type=Path, default=Path(__file__).parents[1] / "results/simulation/simulation_cell_scores.tex")
    args = parser.parse_args()
    rows = read_scores(args.scores)
    if not args.latex.is_file():
        raise FileNotFoundError(args.latex)
    print(f"validated {len(rows)} rows and {len(rows) * len(METRICS)} scores")
    print(f"latex fragment: {args.latex}")


if __name__ == "__main__":
    main()
