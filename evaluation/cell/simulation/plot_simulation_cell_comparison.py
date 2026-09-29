#!/usr/bin/env python3
"""Render source-matched simulation Cell and Gene summary panels."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]

SCENARIOS = ("murk", "branching", "all")
CATEGORIES = (
    ("normal + MURK", "murk"),
    ("normal + branching", "branching"),
    ("hybrid: MURK", "all"),
    ("hybrid: branching", "all"),
)
PLOT_METRICS = (
    ("cbdir", "CBDir", "#d9772a", "D"),
    ("batc", "BATC", "#3a9d5d", "s"),
)
CROSS_LEVEL_MARKER_AREA = 10


def gene_category(row: pd.Series) -> str | None:
    if row.dataset == "murk":
        return "normal + MURK"
    if row.dataset == "branching":
        return "normal + branching"
    if row.dataset == "all":
        if row.gene.startswith(("OneSidedBoost_", "OneSidedBoost2_")):
            return "hybrid: MURK"
        if row.gene.startswith(("MultiForward_", "BoostALL_")):
            return "hybrid: branching"
    return None


def load_summary() -> tuple[pd.DataFrame, pd.DataFrame, list[str], dict]:
    inputs = pd.read_csv(HERE / "inputs.tsv", sep="\t", dtype=str)
    cell = pd.read_csv(HERE / "results/cell_metrics.csv")
    gene = pd.read_csv(REPO / "evaluation/gene/results/gene_cbdir_reference30_all.csv")
    gene = gene[gene.dataset.isin(SCENARIOS) & gene.status.eq("ok")].copy()
    gene["result_path"] = gene.source_h5ad.str.replace(
        "artifact://velocity/", "velocity:", regex=False
    )
    source_map = inputs[["evaluation_dataset", "operational_id", "result_path"]]
    gene = gene.merge(
        source_map, left_on=["dataset", "result_path"],
        right_on=["evaluation_dataset", "result_path"],
        how="left", validate="many_to_one",
    )
    if gene.operational_id.isna().any():
        raise ValueError("a Gene simulation source has no matched Cell input")
    gene["category"] = gene.apply(gene_category, axis=1)
    if gene.category.isna().any():
        raise ValueError("unclassified simulated gene")
    summary = gene.groupby(
        ["category", "dataset", "operational_id"], sort=False
    ).new_cbdir.agg(n_genes="count", gene_cbdir_mean="mean", gene_cbdir_sd="std").reset_index()
    summary = summary.merge(
        cell[["dataset", "method", "cbdir", "batc"]],
        left_on=["dataset", "operational_id"],
        right_on=["dataset", "method"],
        how="left", validate="many_to_one",
    ).drop(columns="method")
    summary = summary.rename(columns={"cbdir": "cell_cbdir", "batc": "cell_batc"})
    expected = {"normal + MURK": 10, "normal + branching": 20,
                "hybrid: MURK": 20, "hybrid: branching": 20}
    if len(summary) != 36 or summary[["cell_cbdir", "cell_batc"]].isna().any().any():
        raise ValueError("expected 36 complete cross-level rows")
    for category, count in expected.items():
        if not summary.loc[summary.category.eq(category), "n_genes"].eq(count).all():
            raise ValueError(f"{category}: unexpected gene count")
    with (REPO / "methods/configurations.tsv").open(encoding="utf-8", newline="") as handle:
        registry = list(csv.DictReader(handle, delimiter="\t"))
    display = {row["operational_id"]: row["manuscript_name"] for row in registry}
    method_order = [row["operational_id"] for row in registry if row["operational_id"] in set(inputs.operational_id)]
    if len(method_order) != 9 or set(method_order) != set(inputs.operational_id):
        raise ValueError("simulation methods differ from the registry subset")
    if len(cell) != 27 or cell.duplicated(["dataset", "method"]).any():
        raise ValueError("expected 27 unique Cell rows")
    summary["display_method"] = summary.operational_id.map(display)
    return cell, summary, method_order, display


def style(ax, show_labels: bool, method_order: list[str], display: dict) -> None:
    ys = np.arange(len(method_order))[::-1]
    for y in ys:
        ax.axhline(y, color="#E9EEF1", lw=0.55, zorder=0)
    ax.axvline(0, color="#9DA7AE", lw=0.8, ls=(0, (2.4, 2.4)), zorder=0)
    ax.set_xlim(-1, 1)
    ax.set_ylim(-0.55, 8.55)
    ax.set_xticks([-1, 0, 1], ["−1", "0", "1"])
    ax.set_yticks(ys, [display[m] for m in method_order])
    ax.tick_params(axis="x", labelsize=7, length=2, colors="#37434B")
    ax.tick_params(axis="y", labelsize=7, length=0, pad=3, colors="#27353E")
    if not show_labels:
        ax.tick_params(axis="y", labelleft=False)
    ax.spines[["left", "right", "top"]].set_visible(False)
    ax.spines["bottom"].set_color("#86929A")


def save(fig, stem: Path) -> None:
    for ext in ("pdf", "svg", "png"):
        fig.savefig(stem.with_suffix(f".{ext}"), dpi=300, facecolor="white")
    svg = stem.with_suffix(".svg")
    svg.write_text(
        "\n".join(line.rstrip() for line in svg.read_text(encoding="utf-8").splitlines()) + "\n",
        encoding="utf-8",
    )
    plt.close(fig)


def plot_cell(cell: pd.DataFrame, method_order: list[str], display: dict, out: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(7.6, 2.8), sharey=True)
    fig.subplots_adjust(left=0.19, right=0.99, top=0.74, bottom=0.23, wspace=0.13)
    titles = {"murk": "Normal + MURK", "branching": "Normal + branching", "all": "Hybrid"}
    for i, dataset in enumerate(SCENARIOS):
        ax = axes[i]
        style(ax, i == 0, method_order, display)
        subset = cell[cell.dataset.eq(dataset)].set_index("method")
        for j, method in enumerate(method_order):
            y = 8 - j
            for column, _, color, marker in PLOT_METRICS:
                ax.scatter(subset.loc[method, column], y + (0.10 if column == "cbdir" else -0.10),
                           color=color, marker=marker, s=23, edgecolor="white", lw=0.45, zorder=3)
        ax.set_title(titles[dataset], fontsize=8, weight="bold", pad=7)
    handles = [plt.Line2D([], [], marker=marker, color=color, lw=0, markersize=5,
                          label=label) for _, label, color, marker in PLOT_METRICS]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.19, 0.985),
               ncol=2, frameon=False, fontsize=7)
    fig.text(0.59, 0.055, "Cell-level score", ha="center", fontsize=7)
    save(fig, out / "simulation_cell_cbdir_batc")


def plot_cross_level(summary: pd.DataFrame, method_order: list[str],
                     display: dict, out: Path) -> None:
    fig, axes = plt.subplots(1, 4, figsize=(9.0, 2.8), sharey=True)
    fig.subplots_adjust(left=0.16, right=0.992, top=0.74, bottom=0.23, wspace=0.12)
    for i, (category, _) in enumerate(CATEGORIES):
        ax = axes[i]
        style(ax, i == 0, method_order, display)
        subset = summary[summary.category.eq(category)].set_index("operational_id")
        for j, method in enumerate(method_order):
            row = subset.loc[method]
            y = 8 - j
            ax.plot([row.gene_cbdir_mean - row.gene_cbdir_sd,
                     row.gene_cbdir_mean + row.gene_cbdir_sd], [y] * 2,
                    color="#3b6ea8", lw=1.3, alpha=0.70, zorder=2)
            ax.scatter(row.cell_batc, y, color="#3a9d5d", marker="s",
                       s=CROSS_LEVEL_MARKER_AREA, edgecolor="white", lw=0.2, zorder=3)
            ax.scatter(row.cell_cbdir, y, color="#d9772a", marker="D",
                       s=CROSS_LEVEL_MARKER_AREA, edgecolor="white", lw=0.2, zorder=4)
            ax.scatter(row.gene_cbdir_mean, y, color="#3b6ea8", marker="o",
                       s=CROSS_LEVEL_MARKER_AREA, edgecolor="white", lw=0.2, zorder=5)
        ax.set_title(category.replace(": ", ":\n"), fontsize=8, weight="bold", pad=7)
    handles = [
        plt.Line2D([], [], marker="o", color="#3b6ea8", lw=1.3, markersize=3.2,
                   label=r"Mean CBDir$_{gene}$ ± gene SD"),
        plt.Line2D([], [], marker="D", color="#d9772a", lw=0, markersize=3.2,
                   label="CBDir"),
        plt.Line2D([], [], marker="s", color="#3a9d5d", lw=0, markersize=3.2,
                   label="BATC"),
    ]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.16, 0.985),
               ncol=3, frameon=False, fontsize=7)
    fig.text(0.58, 0.055, "Directionality score", ha="center", fontsize=7)
    save(fig, out / "simulation_cross_level_cbdir")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=HERE / "figures")
    parser.add_argument("--cross-level-only", action="store_true")
    args = parser.parse_args()
    plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42,
                         "ps.fonttype": 42, "svg.fonttype": "none"})
    cell, summary, method_order, display = load_summary()
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    if not args.cross_level_only:
        summary.to_csv(HERE / "results/cross_level_summary.csv", index=False,
                       float_format="%.17g")
        plot_cell(cell, method_order, display, out)
    plot_cross_level(summary, method_order, display, out)
    print("Wrote cross-level figure" + ("" if args.cross_level_only else " and cell overview"))


if __name__ == "__main__":
    main()
