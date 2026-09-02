#!/usr/bin/env python3
"""Check Dynamo/TSvelo H5AD reconstruction against archived CSV baselines.

The baselines are intentionally test-only. metrics.py reads H5AD objects and
Table S1 during normal evaluation; this script detects an accidental change to
the special H5AD extraction or preprocessing rules.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


MODULE_DIR = Path(__file__).resolve().parent
DEFAULT_CONFIG = MODULE_DIR / "data" / "case3.json"


def read_config(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def read_matrix(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, index_col=0).apply(pd.to_numeric, errors="coerce")


def compare_matrix(name: str, current: Path, historical: Path, tolerance: float) -> None:
    observed = read_matrix(current)
    expected = read_matrix(historical)
    observed.columns = observed.columns.astype(str)
    expected.columns = expected.columns.astype(str)
    if not observed.index.equals(expected.index):
        raise AssertionError(f"{name}: gene order differs from the historical CSV")
    if not observed.columns.equals(expected.columns):
        raise AssertionError(f"{name}: timepoint columns differ from the historical CSV")
    difference = np.nanmax(
        np.abs(observed.to_numpy(dtype=float) - expected.to_numpy(dtype=float))
    )
    if difference > tolerance:
        raise AssertionError(
            f"{name}: max absolute difference {difference:.3e} exceeds {tolerance:.1e}"
        )
    print(f"{name}: pass (max absolute difference {difference:.3e})")


def compare_metrics(
    method: str, metrics_dir: Path, historical: Path, tolerance: float
) -> None:
    observed = pd.read_csv(metrics_dir / "gene_metrics.csv")
    observed = observed[
        observed["Method"].eq(method)
        & observed["Dataset"].eq("rpe1")
        & observed["Parameter"].eq("alpha")
    ].set_index("Gene")
    expected = pd.read_csv(historical).set_index("Gene")
    observed = observed.loc[expected.index]
    if not observed.index.equals(expected.index):
        raise AssertionError(f"{method} metrics: gene order differs")
    if not (
        observed["Valid_Samples"].to_numpy() == expected["Valid_Samples"].to_numpy()
    ).all():
        raise AssertionError(f"{method} metrics: valid-sample counts differ")
    difference = np.nanmax(
        np.abs(
            observed["Spearman_Corr"].to_numpy(dtype=float)
            - expected["Spearman_Corr"].to_numpy(dtype=float)
        )
    )
    if difference > tolerance:
        raise AssertionError(
            f"{method} metrics: max absolute difference {difference:.3e} "
            f"exceeds {tolerance:.1e}"
        )
    print(f"{method} metrics: pass (max absolute difference {difference:.3e})")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--working-root", type=Path, required=True)
    parser.add_argument("--metrics", type=Path, required=True)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--matrix-tolerance", type=float, default=1e-10)
    parser.add_argument("--metric-tolerance", type=float, default=1e-5)
    args = parser.parse_args()

    working_root = args.working_root.resolve()
    metrics_dir = args.metrics.resolve()
    config = read_config(args.config)
    for method, baseline in config["historical_regression"].items():
        compare_matrix(
            f"{method} raw",
            metrics_dir / "matrices" / "rpe1" / f"{method}_alpha_raw.csv",
            working_root / baseline["raw"],
            args.matrix_tolerance,
        )
        compare_matrix(
            f"{method} processed",
            metrics_dir / "trajectories" / "rpe1" / f"{method}_alpha.csv",
            working_root / baseline["processed"],
            args.matrix_tolerance,
        )
        compare_metrics(
            method,
            metrics_dir,
            working_root / baseline["metrics"],
            args.metric_tolerance,
        )


if __name__ == "__main__":
    main()
