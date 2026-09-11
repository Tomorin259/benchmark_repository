#!/usr/bin/env python3
"""Validate the C16 DeepVelo role fields embedded in Cell inputs.tsv."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


MODULE_DIR = Path(__file__).resolve().parent
DEFAULT_INPUTS = MODULE_DIR / "inputs.tsv"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def check_cell_inputs(path: Path) -> int:
    rows = read_tsv(path)
    required = {"operational_id", "evaluation_role", "pred_unspliced", "result_path"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"{path}: missing merged Cell role columns")
    if any(row.get("evaluation_role") != "cell" for row in rows):
        raise ValueError(f"{path}: every Cell input must declare evaluation_role=cell")
    deepvelo = [row for row in rows if row.get("operational_id") == "deepvelo"]
    if len(deepvelo) != 14:
        raise ValueError(f"{path}: expected 14 DeepVelo Cell rows, found {len(deepvelo)}")
    if any(row.get("pred_unspliced") != "false" for row in deepvelo):
        raise ValueError(f"{path}: DeepVelo Cell rows must declare pred_unspliced=false")
    if any("results_with_vu" in row.get("result_path", "") for row in deepvelo):
        raise ValueError(f"{path}: DeepVelo Cell rows still point to the with-vu family")
    if any(
        row.get("pred_unspliced", "")
        for row in rows
        if row.get("operational_id") != "deepvelo"
    ):
        raise ValueError(f"{path}: pred_unspliced is set for a non-DeepVelo row")
    return len(deepvelo)


def check_run_config(path: Path) -> None:
    data = json.loads(path.read_text(encoding="utf-8"))
    expected = {
        "config_id": "C16",
        "operational_id": "deepvelo",
        "evaluation_role": "cell",
        "pred_unspliced": False,
    }
    for key, value in expected.items():
        if data.get(key) != value:
            raise ValueError(f"{path}: {key} does not match the Cell DeepVelo contract")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, default=DEFAULT_INPUTS)
    parser.add_argument("--run-config", type=Path)
    args = parser.parse_args()
    count = check_cell_inputs(args.inputs)
    if args.run_config:
        check_run_config(args.run_config)
    print(f"C16 Cell role contract OK: {count} DeepVelo rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
