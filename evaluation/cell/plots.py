#!/usr/bin/env python3
"""Render the five cell-level metric distributions from compact CSV tables."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


MODULE_DIR = Path(__file__).resolve().parent
DEFAULT_METRICS = MODULE_DIR / "results/all_metrics.csv"
DEFAULT_TRANSPROB = MODULE_DIR / "results/transprob.csv"

METHOD_DISPLAY = {
    "scvelo_deterministic": "scVelo_deterministic",
    "scvelo_stochastic": "scVelo_stochastic",
    "scvelo_dynamic": "scVelo_dynamic",
    "multivelo": "MultiVelo (RNA_only)",
    "unitvelo": "UniTVelo",
    "unitvelo_mode2": "UniTVelo(mode2)",
    "dynamo": "Dynamo(RNA_seq)",
    "TIvelo": "TIvelo",
    "cell2fate": "cell2fate",
    "TFvelo": "TFvelo",
    "veloae": "veloAE",
    "veloVI": "veloVI",
    "DeepKINET": "DeepKINET",
    "TSvelo": "TSvelo",
    "celldancer": "cellDancer",
    "deepvelo": "DeepVelo",
    "SDEvelo": "SDEvelo",
    "deepvelo2022": "DeepVelo_scv",
    "graphvelo_scv": "GraphVelo_scv",
    "graphvelo_dynamo": "GraphVelo_dyn",
}

METHOD_ORDER = [
    "scVelo_deterministic",
    "scVelo_stochastic",
    "scVelo_dynamic",
    "MultiVelo (RNA_only)",
    "UniTVelo",
    "UniTVelo(mode2)",
    "Dynamo(RNA_seq)",
    "TIvelo",
    "cell2fate",
    "TFvelo",
    "veloAE",
    "veloVI",
    "DeepKINET",
    "TSvelo",
    "cellDancer",
    "DeepVelo",
    "SDEvelo",
    "DeepVelo_scv",
    "GraphVelo_scv",
    "GraphVelo_dyn",
]

COLORS = {
    "scVelo_deterministic": "#1F77B4",
    "scVelo_stochastic": "#FF7F0E",
    "Dynamo(RNA_seq)": "#2CA02C",
    "veloAE": "#D62728",
    "scVelo_dynamic": "#9467BD",
    "MultiVelo (RNA_only)": "#8C564B",
    "UniTVelo": "#E377C2",
    "UniTVelo(mode2)": "#7F7F7F",
    "veloVI": "#BCBD22",
    "cell2fate": "#17BECF",
    "SDEvelo": "#0B559F",
    "DeepKINET": "#FA8F12",
    "cellDancer": "#5BA658",
    "DeepVelo": "#C9383B",
    "TFvelo": "#806BA8",
    "TSvelo": "#A05D56",
    "TIvelo": "#D16FA3",
    "DeepVelo_scv": "#5F6B6D",
    "GraphVelo_scv": "#B5BD61",
    "GraphVelo_dyn": "#2E91E5",
}

METRIC_SPECS = [
    ("consistency", "Consistency", "all_metrics", (0.1, 1.05)),
    ("iccoh", "ICCoh", "all_metrics", (0.2, 1.05)),
    ("cbdir", "CBDir", "all_metrics", (-1.05, 1.05)),
    ("transprob", "TransCosine", "transprob", (-0.03, 0.72)),
    ("batc", "BATC", "all_metrics", (-1.05, 1.05)),
]


def load_metric_tables(metrics_path: Path, transprob_path: Path) -> dict[str, pd.DataFrame]:
    all_metrics = pd.read_csv(metrics_path)
    transprob = pd.read_csv(transprob_path)
    for frame, column in ((all_metrics, "all_metrics"), (transprob, "transprob")):
        required = {"dataset", "method", column if column == "transprob" else "consistency"}
        missing = required.difference(frame.columns)
        if missing:
            raise ValueError(f"{column}: missing columns {sorted(missing)}")
    all_metrics["method"] = all_metrics["method"].map(METHOD_DISPLAY).fillna(all_metrics["method"])
    transprob["method"] = transprob["method"].map(METHOD_DISPLAY).fillna(transprob["method"])
    return {"all_metrics": all_metrics, "transprob": transprob}


def ordered_methods(frame: pd.DataFrame) -> list[str]:
    present = set(frame["method"].dropna())
    return [method for method in METHOD_ORDER if method in present]


def draw_metric(axis, frame: pd.DataFrame, column: str, title: str, xlim: tuple[float, float]) -> None:
    methods = ordered_methods(frame)
    grouped = [frame.loc[frame["method"].eq(method), "score"].dropna().to_numpy() for method in methods]
    positions = np.arange(1, len(methods) + 1)
    boxes = axis.boxplot(
        grouped,
        vert=False,
        positions=positions,
        widths=0.55,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "black", "linewidth": 1.2},
        whiskerprops={"linewidth": 1.0},
        capprops={"linewidth": 1.0},
    )
    for patch, method in zip(boxes["boxes"], methods):
        patch.set_facecolor(COLORS.get(method, "#999999"))
        patch.set_alpha(0.65)
    rng = np.random.default_rng(123)
    for position, method, values in zip(positions, methods, grouped):
        if len(values):
            jitter = rng.uniform(-0.13, 0.13, size=len(values))
            axis.scatter(
                values,
                position + jitter,
                s=18,
                color=COLORS.get(method, "#666666"),
                edgecolor="black",
                linewidth=0.35,
                alpha=0.75,
                zorder=3,
            )
            axis.scatter(
                np.mean(values),
                position,
                s=75,
                marker="D",
                color="#FFD700",
                edgecolor="black",
                linewidth=0.9,
                zorder=4,
            )
    axis.axvline(0, color="0.5", linewidth=1.0, alpha=0.65, zorder=0)
    axis.set_yticks(positions, methods)
    axis.tick_params(axis="y", labelsize=12)
    axis.tick_params(axis="x", labelsize=11)
    axis.set_xlim(*xlim)
    axis.set_title(title, fontsize=16, fontweight="bold")
    axis.set_xlabel("Score", fontsize=12)
    axis.grid(axis="x", linestyle="--", alpha=0.35)
    axis.grid(axis="y", visible=False)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)


def make_frame(tables: dict[str, pd.DataFrame], column: str, source: str) -> pd.DataFrame:
    frame = tables[source][["dataset", "method", column]].rename(columns={column: "score"}).copy()
    return frame.dropna(subset=["score"])


def render(metrics_path: Path, transprob_path: Path, output: Path, dpi: int) -> None:
    tables = load_metric_tables(metrics_path, transprob_path)
    output.mkdir(parents=True, exist_ok=True)
    for key, title, source, xlim in METRIC_SPECS:
        figure, axis = plt.subplots(figsize=(4.0, 9.0))
        draw_metric(axis, make_frame(tables, key, source), "score", title, xlim)
        figure.tight_layout()
        for extension in ("pdf", "svg", "png"):
            figure.savefig(output / f"fig2_cell_metric_{key}.{extension}", dpi=dpi)
        plt.close(figure)

    figure, axes = plt.subplots(1, len(METRIC_SPECS), figsize=(15.8, 14.0), sharey=False)
    for axis, (key, title, source, xlim) in zip(axes, METRIC_SPECS):
        draw_metric(axis, make_frame(tables, key, source), "score", title, xlim)
    figure.subplots_adjust(left=0.08, right=0.995, top=0.96, bottom=0.05, wspace=0.75)
    for extension in ("pdf", "svg", "png"):
        figure.savefig(output / f"fig2_cell_metrics_5panel.{extension}", dpi=dpi)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metrics", type=Path, default=DEFAULT_METRICS)
    parser.add_argument("--transprob", type=Path, default=DEFAULT_TRANSPROB)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--dpi", type=int, default=600)
    args = parser.parse_args()
    render(args.metrics.resolve(), args.transprob.resolve(), args.output.resolve(), args.dpi)
    print(f"Wrote cell-level figures to {args.output.resolve()}")


if __name__ == "__main__":
    main()
