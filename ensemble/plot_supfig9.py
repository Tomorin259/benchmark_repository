#!/usr/bin/env python3
"""Render the two reproducible CenterWeighted-4 SupFig9 metric panels."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from evaluation.cell.plots import BENCHMARK_COLORS, METHOD_NAME_MAPPING

SEED = 123
ENSEMBLE_ID = "signed_t_centerweighted4"
ENSEMBLE_LABEL = "CenterWeighted-4"
CONSTITUENTS = {"TIvelo", "veloVI", "deepvelo", "unitvelo"}
DATASET_MAPPING = {"dentate_lamanno": "DentateGyrus_lamanno"}
BASE_DATASET_MAPPING = {"DentateGyrus_lamanno": "dentate_lamanno"}
METRICS = {"cbdir": "CBDir", "batc": "BATC"}
PANEL_STEMS = {"cbdir": "A_centerweighted4_cbdir", "batc": "B_centerweighted4_batc"}


def _check_unique(frame: pd.DataFrame, keys: list[str], label: str) -> None:
    duplicated = frame.duplicated(keys, keep=False)
    if duplicated.any():
        raise ValueError(f"{label} has duplicate keys {keys}")


def load_values(base_csv: Path, ensemble_csv: Path) -> pd.DataFrame:
    base = pd.read_csv(base_csv)
    ensemble = pd.read_csv(ensemble_csv)
    required_base = {"dataset", "method", *METRICS}
    required_ensemble = {"dataset", "variant", *METRICS}
    if missing := required_base - set(base.columns):
        raise ValueError(f"base CSV missing columns: {sorted(missing)}")
    if missing := required_ensemble - set(ensemble.columns):
        raise ValueError(f"ensemble CSV missing columns: {sorted(missing)}")

    base = base[["dataset", "method", *METRICS]].copy()
    _check_unique(base, ["dataset", "method"], "base CSV")
    if base["dataset"].nunique() != 14 or base["method"].nunique() != 20:
        raise ValueError("base CSV must contain 14 datasets and 20 methods")
    if base.groupby("dataset")["method"].nunique().ne(20).any():
        raise ValueError("base CSV must contain all 20 methods for every dataset")
    base["method_label"] = base["method"].map(METHOD_NAME_MAPPING).fillna(base["method"])
    base["dataset_key"] = base["dataset"].replace(BASE_DATASET_MAPPING)
    base["dataset"] = base["dataset"].replace(DATASET_MAPPING)
    _check_unique(base, ["dataset", "method"], "normalized base CSV")

    ensemble = ensemble[ensemble["variant"] == ENSEMBLE_ID].copy()
    if len(ensemble) != 14:
        raise ValueError("ensemble CSV must contain exactly 14 CenterWeighted-4 rows")
    _check_unique(ensemble, ["dataset", "variant"], "ensemble CSV")
    if set(ensemble["dataset"]) != set(base["dataset_key"]):
        raise ValueError("base and CenterWeighted-4 dataset sets differ")
    ensemble["dataset_key"] = ensemble["dataset"]
    ensemble["dataset"] = ensemble["dataset"].replace(DATASET_MAPPING)
    ensemble["method"] = ENSEMBLE_ID
    ensemble["method_label"] = ENSEMBLE_LABEL

    rows = []
    for metric in METRICS:
        current = base[["dataset", "dataset_key", "method", "method_label", metric]].rename(columns={metric: "score"})
        current["metric"] = metric
        extra = ensemble[["dataset", "dataset_key", "method", "method_label", metric]].rename(columns={metric: "score"})
        extra["metric"] = metric
        rows.append(pd.concat([current, extra], ignore_index=True))
    return pd.concat(rows, ignore_index=True)


def common_order(values: pd.DataFrame) -> list[str]:
    base = values[values["method"] != ENSEMBLE_ID]
    metric_ranks = []
    for metric in METRICS:
        means = base[base["metric"] == metric].groupby("method")["score"].mean()
        metric_ranks.append(means.rank(method="average", ascending=False))
    average_rank = pd.concat(metric_ranks, axis=1).mean(axis=1).sort_values(kind="stable")
    return [ENSEMBLE_ID, *average_rank.index.tolist()]


def plot_metric(values: pd.DataFrame, metric: str, order: list[str], output_dir: Path) -> None:
    display_names = {}
    colors = {}
    for method in order:
        if method == ENSEMBLE_ID:
            display_names[method] = ENSEMBLE_LABEL
            colors[ENSEMBLE_LABEL] = "#CC6677"
            continue
        label = values.loc[values["method"] == method, "method_label"].iloc[0]
        if method in CONSTITUENTS:
            label = f"{label} *"
        display_names[method] = label
        source_label = values.loc[values["method"] == method, "method_label"].iloc[0]
        colors[label] = BENCHMARK_COLORS.get(source_label, "#777777")

    frame = values[values["metric"] == metric].copy()
    frame["plot_method"] = frame["method"].map(display_names)
    plot_order = [display_names[item] for item in order]
    means = frame.groupby("plot_method")["score"].mean()

    sns.set_theme(style="ticks", rc={
        "axes.grid": True, "grid.linestyle": "--", "grid.alpha": 0.5,
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    })
    plt.rcParams["svg.fonttype"] = "none"
    plt.rcParams["svg.hashsalt"] = "supfig09-centerweighted4"
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    fig, ax = plt.subplots(figsize=(11, 15.2))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    sns.boxplot(y="plot_method", x="score", data=frame, hue="plot_method", palette=colors,
                order=plot_order, dodge=False, legend=False, width=0.4, linewidth=1.5,
                showfliers=False, boxprops={"alpha": 0.6}, ax=ax)
    np.random.seed(SEED)
    sns.stripplot(y="plot_method", x="score", data=frame, hue="plot_method", palette=colors,
                  order=plot_order, size=5, jitter=0.2, alpha=0.9, linewidth=0.8,
                  edgecolor="black", legend=False, ax=ax)
    for row_id, label in enumerate(plot_order):
        ax.scatter(means[label], row_id, color="#FFD700", s=100, marker="D",
                   edgecolor="black", linewidth=1.2, zorder=10)
    ax.axvline(0, color="gray", linestyle="-", linewidth=2, zorder=0, alpha=0.6)
    ax.set_title(METRICS[metric], fontsize=18, fontweight="bold", pad=10)
    ax.set_xlabel("Score", fontsize=14)
    ax.set_ylabel("")
    ax.tick_params(axis="y", labelsize=14, length=0)
    ax.set_xlim(-1.05, 1.05)
    ax.set_xticks(np.linspace(-1.0, 1.0, 5))
    sns.despine(ax=ax, trim=True, left=True, offset=5)
    for tick in ax.get_yticklabels():
        if tick.get_text() == ENSEMBLE_LABEL:
            tick.set_fontweight("bold")
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = output_dir / PANEL_STEMS[metric]
    fig.tight_layout()
    fig.savefig(
        stem.with_suffix(".svg"),
        format="svg",
        dpi=600,
        facecolor="white",
        metadata={"Date": None},
    )
    plt.close(fig)


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-csv", type=Path, default=root / "evaluation/cell/results/cell_metrics.csv")
    parser.add_argument("--ensemble-csv", type=Path, default=root / "ensemble/results/scores.csv")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=root / "figures/assembly/suppfig09",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    values = load_values(args.base_csv, args.ensemble_csv)
    order = common_order(values)
    for metric in METRICS:
        plot_metric(values, metric, order, args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
