#!/usr/bin/env python3
"""Render real-data Gene/Cell comparisons and common-gene CBDir from CSVs.

Figure 5C/H scatter points average dataset means with equal dataset weights.
The Figure 5H boxplot instead pools valid selected genes from the two datasets.
The S9 panel uses only genes shared by all nine methods within each dataset.
No H5AD, metric recalculation, or gene-output-matched Cell scores are needed.
"""

from __future__ import annotations

import argparse
from importlib.metadata import version
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import seaborn as sns

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
from evaluation.gene.plot_gene_assets import (
    COLORS, NAMES, ORDER, SUMMARY_NAMES, sha256,
)
from evaluation.cell.plots import METHOD_NAME_MAPPING

GROUPS = {
    "MURK-five": ["gastrulation", "capo_pre", "bonemarrow", "mousebrain", "pancreas"],
    "Branching-two": ["DentateGyrus_lamanno", "retina"],
}
SCATTER_ORDER_DATASETS = {
    "MURK-five": ["gastrulation"],
    "Branching-two": ["DentateGyrus_lamanno", "retina"],
}
COMMON_COUNTS = {
    "retina": 6, "DentateGyrus_lamanno": 6, "gastrulation": 6,
    "capo_pre": 7, "bonemarrow": 2, "mousebrain": 4, "pancreas": 4,
}
# Keep the retained S9 Cell panels' historical display label in this panel only.
S9_NAMES = {**NAMES, "scVelo_dyn": "scVelo_dynamic"}
METRICS = ["gene_cbdir", "cell_cbdir", "batc"]
METRIC_STYLE = [
    ("#3B6EA8", "o", r"$\mathrm{CBDir}_{\mathrm{gene}}$"),
    ("#D9772A", "D", "CBDir"),
    ("#3A9D5D", "s", "BATC"),
]


def load_and_validate(paths):
    """Check full-precision summary tables against the canonical gene rows."""
    genes, cell, common, macro, membership = [pd.read_csv(path) for path in paths]
    if genes.duplicated(["dataset", "method", "gene"]).any():
        raise ValueError("Duplicate dataset/method/gene keys in canonical Gene scores")
    real = genes[genes.dataset.isin(COMMON_COUNTS)].copy()
    if set(real.method) != set(ORDER) or set(real.dataset) != set(COMMON_COUNTS):
        raise ValueError("Real Gene scores must cover seven datasets and nine methods")
    finite = np.isfinite(real.new_cbdir)
    if not finite.equals(real.status.eq("ok")):
        raise ValueError("Gene finite-score and status=ok masks disagree")
    valid = real[finite].copy()
    cell["method"] = cell.method.map(METHOD_NAME_MAPPING)
    cell["dataset"] = cell.dataset.replace({"capofibro": "capo_pre"})
    if cell.duplicated(["dataset", "method"]).any():
        raise ValueError("Duplicate primary Cell dataset/method keys")
    cells = cell[cell.dataset.isin(COMMON_COUNTS) & cell.method.isin(NAMES.values())]
    if len(cells) != 63 or not np.isfinite(cells[["cbdir", "batc"]]).all().all():
        raise ValueError("Primary Cell comparison needs 63 finite CBDir/BATC pairs")
    if membership.duplicated(["dataset", "gene"]).any():
        raise ValueError("Duplicate common-gene membership keys")
    if membership.groupby("dataset").size().to_dict() != COMMON_COUNTS:
        raise ValueError("Common-gene counts differ from the seven-dataset contract")
    shared = real.merge(membership, on=["dataset", "gene"], validate="many_to_one")
    if len(shared) != sum(COMMON_COUNTS.values()) * len(ORDER) or not np.isfinite(shared.new_cbdir).all():
        raise ValueError("Common genes must have valid scores for all nine methods")
    expected = shared.groupby(["dataset", "method"]).new_cbdir.agg(
        n_common_genes="size", mean="mean", sd="std")
    keys = ["dataset", "method"]
    if len(common) != 63 or common.duplicated(keys).any():
        raise ValueError("Common-gene dataset summary must contain 63 unique rows")
    pd.testing.assert_frame_equal(
        common.set_index(keys)[["n_common_genes", "mean", "sd"]].sort_index(),
        expected.sort_index(), check_exact=False, rtol=0, atol=1e-14)
    calculated = common.groupby("method")["mean"].mean().sort_index()
    if len(macro) != 9 or macro.method.duplicated().any() or not macro.n_datasets.eq(7).all():
        raise ValueError("Common-gene macro summary must contain nine seven-dataset means")
    np.testing.assert_allclose(macro.set_index("method")["mean"].sort_index(),
                               calculated, rtol=0, atol=1e-14)
    ranks = calculated.rank(ascending=False, method="min").astype(int)
    np.testing.assert_array_equal(macro.set_index("method")["rank"].sort_index(), ranks)
    gene_means = valid.groupby(["dataset", "method"]).new_cbdir.agg(
        value="mean", n_valid_genes="size").reset_index()
    gene_means["metric"] = "gene_cbdir"
    primary = cells.rename(columns={"cbdir": "cell_cbdir"}).melt(
        id_vars=["dataset", "method"], value_vars=["cell_cbdir", "batc"],
        var_name="metric", value_name="value")
    primary["method"] = primary.method.map({name: key for key, name in NAMES.items()})
    dataset_means = pd.concat([gene_means, primary], ignore_index=True)
    return real, common, macro, dataset_means


def plot_scatter(frame, order, output, dpi):
    """Use adjacent gene-distribution order and equal dataset metric weights."""
    plt.rcParams.update(plt.rcParamsDefault)
    fig, ax = plt.subplots(figsize=(7.5, 3.15))
    for metric, (color, marker, label) in zip(METRICS, METRIC_STYLE):
        values = frame[frame.metric == metric].set_index("method").loc[order, "value"]
        ax.scatter(values, np.arange(len(order)), c=color, marker=marker,
                   s=24, edgecolors=color, linewidths=0.5, label=label, zorder=3)
    ax.set_yticks(np.arange(len(order)), [SUMMARY_NAMES[method] for method in order])
    ax.set_ylim(len(order) - 0.4, -0.6)
    ax.set_xlim(-0.25, 1)
    ax.set_xticks([-0.25, 0, 0.5, 1])
    ax.set_xlabel("score", fontsize=10)
    ax.tick_params(axis="both", labelsize=9)
    ax.tick_params(axis="y", length=0, pad=6)
    ax.xaxis.grid(True, color="0.87", linewidth=0.6)
    ax.axvline(0, color="0.35", linewidth=0.9, zorder=1)
    ax.legend(loc="upper right", bbox_to_anchor=(1, -0.15), ncol=3,
              frameon=False, fontsize=8, handletextpad=0.3, columnspacing=0.9,
              borderaxespad=0, handlelength=0.9)
    fig.subplots_adjust(left=0.205, right=0.975, bottom=0.27, top=0.89)
    for suffix in ["png", "pdf", "svg"]:
        fig.savefig(output.with_suffix("." + suffix), dpi=dpi,
                    facecolor="white", transparent=False)
    plt.close(fig)


def plot_branching(real, output, dpi):
    """Preserve the established pooled-gene boxplot/heatmap, including NaNs."""
    selected = real[real.dataset.isin(GROUPS["Branching-two"])].copy()
    if len(selected) != 180 or selected.gene.nunique() != 20:
        raise ValueError("The branching panel expects 180 selected gene rows and 20 genes")
    frame = selected[np.isfinite(selected.new_cbdir)].copy()
    frame["Method"] = frame.method.map(SUMMARY_NAMES)
    frame = frame.rename(columns={"gene": "Gene", "new_cbdir": "Overall_Score"})
    if len(frame) != 165:
        raise ValueError("The branching panel expects 165 valid selected-gene scores")
    means = frame.groupby("Method").Overall_Score.mean()
    order = means.sort_values(ascending=False).index.tolist()
    sns.set_theme(style="whitegrid", font_scale=1.1)
    fig, axes = plt.subplots(1, 2, figsize=(33, 6),
                            gridspec_kw={"width_ratios": [2.05, 2.70]})
    sns.boxplot(x="Method", y="Overall_Score", data=frame, order=order,
                palette=COLORS, ax=axes[0])
    for index, method in enumerate(order):
        axes[0].scatter(index, means[method], color="#FFD700", s=48, marker="D",
                        edgecolor="black", linewidth=0.6, zorder=10)
    axes[0].set_title("")
    axes[0].set_ylabel(r"$\mathrm{CBDir}_{\mathrm{gene}}$")
    axes[0].set_xlabel("")
    axes[0].tick_params(axis="x", rotation=0, labelsize=11, pad=10)
    heatmap = frame.pivot(index="Method", columns="Gene", values="Overall_Score").loc[order]
    sns.heatmap(heatmap, mask=heatmap.isna(), annot=True, cmap="RdBu_r", center=0,
                fmt=".2f", ax=axes[1], linewidths=0, linecolor="none",
                annot_kws={"size": 8}, cbar_kws={"label": "", "pad": 0.015})
    axes[1].set_facecolor("white")
    axes[1].grid(False)
    axes[1].legend(handles=[mpatches.Patch(facecolor="white", edgecolor="gray", label="Not Calculated")],
                   loc="lower right", bbox_to_anchor=(1, -0.15), frameon=False)
    axes[1].set_title("")
    axes[1].set_ylabel("")
    axes[1].set_xlabel("")
    axes[1].tick_params(axis="x", rotation=0, labelsize=9, pad=8)
    plt.subplots_adjust(bottom=0.18, wspace=0.18)
    plt.tight_layout(rect=[0, 0.05, 1, 1])
    for suffix in ["png", "pdf", "svg"]:
        fig.savefig(output.with_suffix("." + suffix), dpi=dpi, facecolor="white")
    plt.close(fig)
    selected.to_csv(output.parent / "branching_gene_scores.csv", index=False)
    means.rename("pooled_gene_mean").sort_values(ascending=False).to_csv(
        output.parent / "branching_pooled_means.csv")


def plot_common(common, macro, output, dpi):
    """Seven dataset means per method; diamonds are equal-dataset macro means."""
    order = macro.sort_values("mean", ascending=False).method.tolist()
    sns.set_theme(style="whitegrid", font_scale=1.0)
    fig, ax = plt.subplots(figsize=(3.2, 2.7))
    frame = common.copy()
    frame["Method"] = frame.method.map(S9_NAMES)
    colors = {S9_NAMES[method]: COLORS[SUMMARY_NAMES[method]] for method in ORDER}
    sns.boxplot(data=frame, x="mean", y="Method", order=[S9_NAMES[m] for m in order],
                palette=colors, width=0.42, showfliers=False, linewidth=0.65,
                saturation=0.8, ax=ax, boxprops={"alpha": 0.65})
    rng = np.random.default_rng(123)
    frame = frame.sort_values(["method", "dataset"]).reset_index(drop=True)
    frame["plot_y"] = frame.method.map({method: i for i, method in enumerate(order)})
    frame["plot_y"] += rng.uniform(-0.12, 0.12, len(frame))
    for method in order:
        points = frame[frame.method == method]
        ax.scatter(points["mean"], points.plot_y, s=7,
                   color=COLORS[SUMMARY_NAMES[method]], edgecolor="black",
                   linewidth=0.3, alpha=0.8, zorder=4)
    means = macro.set_index("method").loc[order, "mean"]
    ax.scatter(means, np.arange(len(order)), color="#FFD700", marker="D",
               edgecolor="black", linewidth=0.35, s=16, zorder=5)
    ax.axvline(0, color="0.7", linewidth=0.95, zorder=1)
    ax.set_xlim(-1.05, 1.05)
    ax.set_xticks([-1, -0.5, 0, 0.5, 1])
    ax.set_ylim(len(order) - 0.5, -0.5)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("CBDir-gene", fontsize=10, fontweight="bold", pad=6)
    ax.tick_params(axis="both", labelsize=6.8)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", color="0.15", width=0.8, length=2.5, pad=2.5)
    ax.yaxis.grid(False)
    ax.xaxis.grid(True, linestyle="--", color="0.9", linewidth=0.5)
    sns.despine(ax=ax, left=True, bottom=False)
    ax.spines["bottom"].set_color("0.15")
    ax.spines["bottom"].set_linewidth(0.8)
    fig.subplots_adjust(left=0.39, right=0.98, bottom=0.10, top=0.86)
    for suffix in ["png", "pdf", "svg"]:
        fig.savefig(output.with_suffix("." + suffix), dpi=dpi,
                    facecolor="white", transparent=False)
    plt.close(fig)
    frame.to_csv(output.parent / "common_gene_dataset_points.csv", index=False)
    macro.sort_values("rank").to_csv(output.parent / "common_gene_macro_means.csv", index=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--gene-scores", type=Path, default=HERE / "results/gene_cbdir_reference30_all.csv")
    parser.add_argument("--cell-scores", type=Path, default=ROOT / "evaluation/cell/results/cell_metrics.csv")
    parser.add_argument("--common-summary", type=Path, default=HERE / "results/common_genes/dataset_summary.csv")
    parser.add_argument("--common-macro", type=Path, default=HERE / "results/common_genes/macro_summary.csv")
    parser.add_argument("--common-membership", type=Path, default=HERE / "results/common_genes/membership.csv")
    parser.add_argument("--dpi", type=int, default=600, help="Branching and S9 PNG resolution")
    parser.add_argument("--scatter-dpi", type=int, default=300, help="300 gives the original 2250 x 945 size")
    args = parser.parse_args()
    if args.dpi <= 0 or args.scatter_dpi <= 0:
        parser.error("DPI values must be positive")
    output = args.output_root.expanduser().resolve()
    canonical = [HERE / "figures", ROOT / "figures"]
    if output == ROOT or output == ROOT / "evaluation" or any(
            output == folder or folder in output.parents for folder in canonical):
        parser.error("Render to a separate output directory")
    if output.exists() and any(output.iterdir()):
        parser.error("--output-root must be empty; use a new directory")
    paths = [args.gene_scores, args.cell_scores, args.common_summary,
             args.common_macro, args.common_membership]
    real, common, macro, dataset_means = load_and_validate(paths)
    output.mkdir(parents=True, exist_ok=True)
    plotted, scatter_orders = [], {}
    for group, datasets in GROUPS.items():
        rows = dataset_means[dataset_means.dataset.isin(datasets)].copy()
        counts = rows.groupby(["method", "metric"]).size()
        if len(counts) != len(ORDER) * len(METRICS) or not counts.eq(len(datasets)).all():
            raise ValueError(f"{group}: missing dataset/method/metric entries")
        means = rows.groupby(["method", "metric"]).value.mean().reset_index()
        means.insert(0, "group", group)
        means["n_datasets"] = len(datasets)
        means["weighting"] = "equal_dataset_weight"
        means["display_method"] = means.method.map(SUMMARY_NAMES)
        adjacent = real[real.dataset.isin(SCATTER_ORDER_DATASETS[group]) & np.isfinite(real.new_cbdir)]
        order = adjacent.groupby("method").new_cbdir.mean().sort_values(ascending=False).index.tolist()
        scatter_orders[group] = order
        means["plot_y"] = means.method.map({method: i for i, method in enumerate(order)})
        plotted.append(means)
        name = "murk" if group == "MURK-five" else "branching"
        plot_scatter(means, order, output / f"{name}_metric_scatter_ordered", args.scatter_dpi)
    pd.concat(plotted, ignore_index=True).to_csv(output / "plotted_scores.csv", index=False)
    dataset_means.to_csv(output / "cross_dataset_means.csv", index=False)
    plot_branching(real, output / "branching_combined_reference30_label_order", args.dpi)
    plot_common(common, macro, output / "common_gene_cbdir", args.dpi)
    canonical_paths = [HERE / "results/gene_cbdir_reference30_all.csv",
                       ROOT / "evaluation/cell/results/cell_metrics.csv",
                       HERE / "results/common_genes/dataset_summary.csv",
                       HERE / "results/common_genes/macro_summary.csv",
                       HERE / "results/common_genes/membership.csv"]
    provenance = {
        "renderer": "evaluation/gene/plot_real_comparisons.py",
        "renderer_sha256": sha256(__file__),
        "source_status": "ACTUAL_RUN",
        "source_tables": {canonical.relative_to(ROOT).as_posix(): sha256(actual)
                          for canonical, actual in zip(canonical_paths, paths)},
        "style_sources": {path.relative_to(ROOT).as_posix(): sha256(path) for path in
                          [HERE / "plot_gene_assets.py", ROOT / "evaluation/cell/plots.py",
                           ROOT / "methods/configurations.tsv"]},
        "software_versions": {package: version(package) for package in
                              ["numpy", "pandas", "matplotlib", "seaborn"]},
        "scatter_method_order": scatter_orders,
        "scatter_ordering": "descending pooled-gene mean in adjacent distribution panel; not the scatter metric mean",
        "scatter_order_datasets": SCATTER_ORDER_DATASETS,
        "scatter_axis_limits": [-0.25, 1],
        "scatter_dataset_groups": GROUPS,
        "scatter_weighting": "within-method valid selected-gene means, then equal dataset weights; primary Cell means",
        "cell_identifier_aliases": "primary Cell operational methods use evaluation/cell/plots.py; capofibro dataset joins capo_pre",
        "branching_box_heatmap_weighting": "pooled valid selected genes from two datasets; yellow diamonds pooled means",
        "common_gene_counts": COMMON_COUNTS,
        "common_gene_display_aliases": {"scVelo_dyn": "scVelo_dynamic"},
        "common_gene_weighting": "seven dataset points per method; yellow diamonds equal-dataset means",
        "common_gene_boxplot": "median line, interquartile box, 1.5-IQR whiskers; all seven dataset points shown",
        "jitter": "NumPy default_rng(123), display only",
        "dpi": args.dpi, "scatter_dpi": args.scatter_dpi,
        "outputs": {path.name: sha256(path) for path in sorted(output.iterdir())},
    }
    (output / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(f"Rendered four real-data summary panels into {output}")


if __name__ == "__main__":
    main()
