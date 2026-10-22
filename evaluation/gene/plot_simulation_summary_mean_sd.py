#!/usr/bin/env python3
"""Render gene-level simulation summaries.

Figure 4D is rendered by the Cell simulation module.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image


CATEGORIES = [
    ("normal + MURK", "MURK"),
    ("normal + branching", "branching"),
    ("hybrid: MURK", "Hybrid:\nMURK"),
    ("hybrid: branching", "Hybrid:\nbranching"),
]
METHODS = [
    "scVelo_dyn", "MultiVelo", "UniTVelo", "TIvelo", "veloVI",
    "cellDancer", "DeepVelo", "DeepVelo_scv", "Graphvelo_scv",
]
DISPLAY = {"MultiVelo": "MultiVelo (RNA_only)", "Graphvelo_scv": "GraphVelo_scv"}
METRICS = [
    ("CBDir", r"$\mathrm{CBDir}_{\mathrm{gene}}$", "#246782", "o", -0.12),
    ("CosSim", "CosSim", "#D47B3A", "D", 0.12),
]


def category(row: pd.Series) -> str | None:
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


def summarize(repo: Path) -> pd.DataFrame:
    cbdir = pd.read_csv(repo / "evaluation/gene/results/gene_cbdir_reference30_all.csv")
    cbdir = cbdir[cbdir.dataset.isin(("murk", "branching", "all"))]
    assert set(cbdir.status) == {"ok"}
    compact = pd.read_csv(repo / "evaluation/gene/results/simulation/gene_cbdir_cossim_reference.csv")
    compact = compact.rename(columns={"Dataset": "dataset", "Gene": "gene", "Method": "method"})
    raw = cbdir[["dataset", "method", "gene", "new_cbdir"]].merge(
        compact[["dataset", "method", "gene", "CosSim"]],
        on=["dataset", "method", "gene"], how="left", validate="one_to_one",
    )
    assert raw.CosSim.notna().all() and not raw.duplicated(["dataset", "method", "gene"]).any()
    raw = raw.rename(columns={"new_cbdir": "CBDir"})
    raw["category"] = raw.apply(category, axis=1)
    raw = raw[raw.category.notna() & raw.method.isin(METHODS)]
    long = raw.melt(
        id_vars=["category", "method", "gene"],
        value_vars=["CBDir", "CosSim"], var_name="metric", value_name="score",
    )
    summary = long.groupby(["category", "method", "metric"], sort=False)["score"].agg(
        n="count", mean="mean", sd="std",
    ).reset_index()
    assert len(summary) == 4 * 9 * 2 and summary.sd.notna().all()
    return summary


def plot(summary: pd.DataFrame, out: Path) -> None:
    mpl.rcParams.update({
        "font.family": "DejaVu Sans", "pdf.fonttype": 42,
        "ps.fonttype": 42, "svg.fonttype": "none",
    })
    fig, axes = plt.subplots(1, 4, figsize=(7.5, 2.20), sharey=True)
    fig.subplots_adjust(left=0.18, right=0.987, top=0.73, bottom=0.22, wspace=0.12)
    balanced = summary.groupby("method", as_index=False)["mean"].mean()
    method_order = balanced.sort_values(
        ["mean", "method"], ascending=[False, True]
    )["method"].tolist()
    y_base = np.arange(len(method_order))[::-1]

    for col, (category_key, title) in enumerate(CATEGORIES):
        ax = axes[col]
        ax.axvline(0, color="#A7AEB5", lw=0.85, ls=(0, (2.2, 2.2)), zorder=0)
        for y in y_base:
            ax.axhline(y, color="#ECF0F2", lw=0.55, zorder=0)
        for metric, _, color, marker, offset in METRICS:
            rows = summary[(summary.category == category_key) & (summary.metric == metric)]
            by_method = rows.set_index("method")
            for i, method in enumerate(method_order):
                row = by_method.loc[method]
                y = y_base[i] + offset
                ax.plot([row["mean"] - row.sd, row["mean"] + row.sd], [y, y],
                        color=color, lw=1.45, solid_capstyle="round", alpha=0.75, zorder=2)
                ax.scatter([row["mean"]], [y], s=16, facecolor=color,
                           edgecolor="white", linewidth=0.4, marker=marker, zorder=3)
        ax.set_xlim(-1, 1)
        ax.set_ylim(-0.55, 8.55)
        ax.set_xticks([-1, 0, 1], labels=["−1", "0", "1"])
        ax.tick_params(axis="x", labelsize=6.1, length=2.1, pad=2.5, colors="#353B40")
        ax.tick_params(axis="y", length=0, pad=3)
        ax.set_title(title, fontsize=7.1, weight="bold", color="#26333B", pad=7)
        ax.spines[["top", "left", "right"]].set_visible(False)
        ax.spines["bottom"].set_color("#7A848B")
        ax.spines["bottom"].set_linewidth(0.6)
        if col == 0:
            ax.set_yticks(y_base, [DISPLAY.get(method, method) for method in method_order])
            ax.tick_params(axis="y", labelleft=True, labelsize=5.9, colors="#26333B")
        else:
            ax.tick_params(axis="y", labelleft=False)

    handles = [
        mpl.lines.Line2D([0], [0], marker=marker, color=color, linewidth=1.45,
                         markersize=3.8, markerfacecolor=color,
                         markeredgecolor="white", label=label)
        for _, label, color, marker, _ in METRICS
    ]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.18, 0.98),
               ncol=2, frameon=False, fontsize=6.2, handlelength=1.5,
               columnspacing=1.25, handletextpad=0.4)
    fig.text(0.59, 0.055, "Gene-level score", fontsize=6.3,
             color="#26333B", ha="center")
    png = out / "summary_dot_interval.png"
    fig.savefig(png, dpi=300, facecolor="white")
    fig.savefig(out / "summary_dot_interval.svg", facecolor="white")
    plt.close(fig)

    with Image.open(png) as source:
        assert source.size == (2250, 660)
        canvas = Image.new("RGB", (2250, 744), "white")
        canvas.paste(source.convert("RGB"), (0, 42))
        canvas.save(out / "summary_dot_interval_ppt.png", dpi=(300, 300))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    repo = args.repo_root.resolve()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    summary = summarize(repo)
    summary.to_csv(out / "summary_values.csv", index=False, float_format="%.17g")
    plot(summary, out)


if __name__ == "__main__":
    main()
