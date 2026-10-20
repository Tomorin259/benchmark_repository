#!/usr/bin/env python3
"""Assemble Ncald and FAM111A with the accepted Fech 6+3 geometry."""
from __future__ import annotations

from pathlib import Path
import argparse
import json
import sys
import hashlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np
import pandas as pd
from PIL import Image


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
from evaluation.gene.plot_gene_assets import load_scores, read_phase_inputs, ORDER, plot_phase_portrait
METHODS = [
    "scVelo_dyn", "MultiVelo", "UniTVelo", "TIvelo", "veloVI",
    "cellDancer", "DeepVelo", "DeepVelo_scv", "Graphvelo_scv",
]
DISPLAY_LABELS = {"Graphvelo_scv": "GraphVelo_scv"}
CASES = {
    "Fech": {"source": "Erythroid2", "target": "Erythroid3", "clusters": [
        ("#2ca02c", "Erythroid1"), ("#d62728", "Erythroid2"), ("#9467bd", "Erythroid3")
    ]},
    "Ncald": {"source": "ImmGranule1", "target": "ImmGranule2", "clusters": [
        ("#e6afb9", "Nbl2"), ("#4a6fe3", "ImmGranule1"),
        ("#8595e1", "ImmGranule2"), ("#bb7784", "Granule")
    ]},
    "FAM111A": {"source": "G1", "target": "S", "clusters": [
        ("#1f77b4", "G1"), ("#2ca02c", "S"), ("#ff7f0e", "G2M")
    ]},
}
FIG_W, FIG_H = 7.5, 2.55



def metric_rows(source_data: pd.DataFrame, gene: str) -> pd.DataFrame:
    q = source_data[(source_data.artifact_type == "scatter") & (source_data.gene == gene)].copy()
    if sorted(q.method.unique()) != sorted(METHODS) or len(q) != len(METHODS):
        raise ValueError(f"{gene}: expected one scatter row per method")
    q["gene_cbdir"] = pd.to_numeric(q["gene_cbdir"], errors="raise")
    q["cell"] = pd.to_numeric(q["cell_cbdir"], errors="raise")
    q["batc"] = pd.to_numeric(q["selected_edge_batc_candidate_mean"], errors="raise")
    return q.set_index("method").loc[METHODS].reset_index()


def draw_metric(ax, data: pd.DataFrame):
    y = list(range(len(data)))[::-1]
    ax.axvline(0, color="#555555", lw=.65, zorder=1)
    for yi, (_, row) in zip(y, data.iterrows()):
        ax.scatter(row.gene_cbdir, yi, s=8, marker="o", color="#3B6EA8", zorder=3)
        ax.scatter(row.cell, yi, s=9, marker="D", color="#D9772A", zorder=3)
        ax.scatter(row.batc, yi, s=9, marker="s", color="#3A9D5D", zorder=3)
    ax.set_xlim(-1, 1); ax.set_ylim(-.6, len(data) - .4)
    ax.set_yticks(y); ax.set_yticklabels([DISPLAY_LABELS.get(method, method) for method in METHODS], fontsize=5.2)
    ax.set_xticks([-1, -.5, 0, .5, 1]); ax.tick_params(axis="x", labelsize=6.2, length=2)
    ax.tick_params(axis="y", length=0, pad=2); ax.grid(axis="x", color="#dddddd", lw=.35)
    ax.set_axisbelow(True); ax.set_xlabel("score", fontsize=7.0, labelpad=2)
    for spine in ax.spines.values(): spine.set_linewidth(.55)


def _transition_arrow(ax, start, end, *, color, selected=False, rad=0.0):
    """Draw a compact transition arrow in axes-fraction coordinates."""
    ax.add_patch(FancyArrowPatch(
        start, end,
        transform=ax.transAxes,
        arrowstyle="-|>",
        mutation_scale=5.2 if selected else 4.6,
        linewidth=1.05 if selected else .75,
        color=color,
        shrinkA=4.5,
        shrinkB=4.5,
        connectionstyle=f"arc3,rad={rad}",
        zorder=3,
    ))


def _transition_node(ax, xy, color, label, *, fontsize=5.0, label_y=None,
                     label_va="top"):
    ax.plot(
        xy[0], xy[1],
        transform=ax.transAxes,
        marker="o",
        linestyle="None",
        markerfacecolor=color,
        markeredgecolor="black",
        markeredgewidth=.55,
        markersize=4.2,
        zorder=4,
    )
    if label_y is None:
        label_y = .055 if xy[1] < .25 else .035
    ax.text(
        xy[0], label_y,
        label,
        transform=ax.transAxes,
        ha="center",
        va=label_va,
        fontsize=fontsize,
        color="#222222",
        zorder=5,
    )


def _selected_box(ax, x0, y0, x1, y1):
    ax.add_patch(FancyBboxPatch(
        (x0, y0), x1 - x0, y1 - y0,
        transform=ax.transAxes,
        boxstyle="round,pad=0.008,rounding_size=0.02",
        fill=False,
        edgecolor="#c2410c",
        linewidth=.9,
        linestyle=(0, (3, 2)),
        zorder=2,
    ))


def draw_transition_graph(ax, case):
    """Draw the full biological relation below the cluster legend.

    The old ``edge:`` key is intentionally not rendered. The dashed box is
    the only indicator of the evaluated transition.
    """
    colors = {label: color for color, label in case["clusters"]}
    source = case["source"]
    target = case["target"]

    if case is CASES["FAM111A"]:
        positions = {"G1": (.31, .27), "S": (.69, .27), "G2M": (.50, .13)}
        _transition_arrow(ax, positions["G1"], positions["S"], color="#c2410c", selected=True)
        _transition_arrow(ax, positions["S"], positions["G2M"], color="#6b7280")
        _transition_arrow(ax, positions["G2M"], positions["G1"], color="#6b7280")
        _selected_box(ax, .27, .205, .73, .335)
        _transition_node(ax, positions["G1"], colors["G1"], "G1",
                         fontsize=5.3, label_y=.35, label_va="bottom")
        _transition_node(ax, positions["S"], colors["S"], "S",
                         fontsize=5.3, label_y=.35, label_va="bottom")
        _transition_node(ax, positions["G2M"], colors["G2M"], "G2M", label_y=.035)
        return

    if case is CASES["Ncald"]:
        labels = ["Nbl2", "ImmGranule1", "ImmGranule2", "Granule"]
        display_labels = {"Nbl2": "Nbl2", "ImmGranule1": "ImmG1",
                          "ImmGranule2": "ImmG2", "Granule": "Granule"}
        positions = dict(zip(labels, ((.14, .22), (.38, .22), (.62, .22), (.86, .22))))
        for left, right in zip(labels[:-1], labels[1:]):
            selected = left == source and right == target
            _transition_arrow(
                ax, positions[left], positions[right],
                color="#c2410c" if selected else "#6b7280",
                selected=selected,
            )
        _selected_box(ax, .33, .16, .67, .285)
        for label in labels:
            _transition_node(ax, positions[label], colors[label], display_labels[label],
                             fontsize=4.6, label_y=.05)
        return

    labels = ["Erythroid1", "Erythroid2", "Erythroid3"]
    positions = dict(zip(labels, ((.22, .22), (.50, .22), (.78, .22))))
    for left, right in zip(labels[:-1], labels[1:]):
        selected = left == source and right == target
        _transition_arrow(
            ax, positions[left], positions[right],
            color="#c2410c" if selected else "#6b7280",
            selected=selected,
        )
    _selected_box(ax, .45, .16, .83, .285)
    for label in labels:
        _transition_node(ax, positions[label], colors[label], label,
                         fontsize=3.6, label_y=.05)


def assemble(gene, source_data, output_dir):
    case = CASES[gene]
    data = metric_rows(source_data, gene)
    fig = plt.figure(figsize=(FIG_W, FIG_H), facecolor="white")
    left, right, gap = .006, .998, .006
    width = (right - left - 5 * gap) / 6
    height = width * FIG_W / FIG_H * (900 / 975)
    bottom = .035
    row_y = [bottom, 1.0 - height - .025]
    phase_dir = output_dir / "phase_panels" / gene
    for i, method in enumerate(METHODS):
        row, col = (1, i) if i < 6 else (0, i - 6)
        x = left + col * (width + gap)
        ax = fig.add_axes([x, row_y[row], width, height])
        shift = 6 if gene == "Ncald" else 0
        ax.imshow(np.asarray(Image.open(phase_dir / f"phase_{method}.png").convert("RGB")), interpolation="lanczos", aspect="equal")
        ax.axis("off")

    legend_box = fig.add_axes([left + 3 * (width + gap), row_y[0], width, height])
    legend_box.axis("off")
    handles = [Line2D([0], [0], marker="o", color="none", markerfacecolor=color,
                      markeredgecolor="black", markeredgewidth=.6,
                      markersize=4.2, label=label)
               for color, label in case["clusters"]]
    legend_box.legend(handles=handles, loc="upper center", bbox_to_anchor=(.5, .80),
                      ncol=1, frameon=False, fontsize=5.1, handletextpad=.25,
                      borderaxespad=0)
    # The former single-edge key is replaced by the complete biological
    # relation; the dashed box marks the transition used for evaluation.
    draw_transition_graph(legend_box, case)

    box = fig.add_axes([left + 4 * (width + gap), row_y[0], 2 * width + gap, height])
    box.axis("off")
    draw_metric(box.inset_axes([.16, .31, .82, .62]), data)
    metric_handles = [
        Line2D([0], [0], marker="o", linestyle="None", color="none",
               markerfacecolor="#3B6EA8", markeredgecolor="#3B6EA8",
               markersize=3.2, label=r"$\mathrm{CBDir}_{\mathrm{gene}}$"),
        Line2D([0], [0], marker="D", linestyle="None", color="none",
               markerfacecolor="#D9772A", markeredgecolor="#D9772A",
               markersize=3.2, label="CBDir"),
        Line2D([0], [0], marker="s", linestyle="None", color="none",
               markerfacecolor="#3A9D5D", markeredgecolor="#3A9D5D",
               markersize=3.2, label="BATC"),
    ]
    fig.legend(handles=metric_handles, loc="lower right", bbox_to_anchor=(.985, .018),
               ncol=3, frameon=False, fontsize=4.8, handletextpad=.22,
               columnspacing=.55, borderaxespad=0, markerscale=1)

    out_dir = output_dir / "selected_edges"
    out_dir.mkdir(exist_ok=True)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(out_dir / f"{gene}_matched_edge.{ext}", dpi=600 if ext != "png" else 300, facecolor="white")
    plt.close(fig)
    data.to_csv(out_dir / f"{gene}_plotted_scores.csv", index=False, float_format="%.17g")
    source_data[source_data.gene == gene].to_csv(out_dir / f"{gene}_source_data.csv", index=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    output = args.output_root.resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError("Use a fresh output directory")
    output.mkdir(parents=True, exist_ok=True)
    scores = load_scores()
    edges = pd.read_csv(HERE / "results/gene_cbdir_edge_details_reference30_all.csv")
    controls = pd.read_csv(HERE / "results/selected_edge_cell_controls.csv")
    datasets = {"Fech": "gastrulation", "Ncald": "DentateGyrus_lamanno", "FAM111A": "capo_pre"}
    input_records = {}
    for gene, case in CASES.items():
        dataset = datasets[gene]
        arrays, inputs = read_phase_inputs(dataset, scores, args.artifact_root.resolve(), [gene])
        input_records[gene] = inputs
        folder = output / "phase_panels" / gene
        folder.mkdir(parents=True)
        source = []
        for method, adata in zip(ORDER, arrays):
            selected = adata[adata.obs["clusters"].astype(str).isin([x[1] for x in case["clusters"]])].copy()
            edge = edges[(edges.dataset == dataset) & (edges.gene == gene) & (edges.method == method)
                         & (edges.source == case["source"]) & (edges.target == case["target"])]
            cell = controls[(controls.gene == gene) & (controls.method == method)]
            if len(edge) != 1 or len(cell) != 1:
                raise ValueError(f"Selected transition is not unique: {gene}/{method}")
            row = {**cell.iloc[0].to_dict(), "gene_cbdir": float(edge.iloc[0].new_edge_score),
                   "artifact_type": "scatter", "dataset": dataset}
            source.append(row)
            fig, ax = plt.subplots(figsize=(3.25, 3.00))
            plot_phase_portrait(selected, gene, ax=ax, color_key="clusters", legend_loc=None,
                                s=12, alpha=0.6, arrow_grid=(25, 25), arrow_density_percentile=15.0,
                                arrow_length_mode="log1p", arrow_reference_quantile=0.90,
                                arrow_reference_length=0.05, arrow_max_length=0.10)
            ax.set_title(""); ax.set_xlabel(""); ax.set_ylabel(""); ax.axis("off")
            ax.set_position([0.12, 0.015, 0.76, 0.82])
            for y, text in [(0.994, DISPLAY_LABELS.get(method, method)),
                            (0.943, r"$\mathrm{CBDir}_{\mathrm{gene}}$: " + f"{row['gene_cbdir']:.4f}"),
                            (0.878, f"CBDir: {float(row['cell_cbdir']):.4f}")]:
                fig.text(0.5, y, text, ha="center", va="top", fontsize=10.5)
            fig.savefig(folder / f"phase_{method}.png", dpi=300, facecolor="white")
            plt.close(fig)
        assemble(gene, pd.DataFrame(source), output)
    paths = [HERE / "results/gene_cbdir_edge_details_reference30_all.csv", HERE / "results/selected_edge_cell_controls.csv", Path(__file__)]
    (output / "provenance.json").write_text(json.dumps({
        "inputs": input_records, "sources": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        "display": "native method coordinates; selected-cluster phase portraits; original six-plus-three layout",
        "scoring": "current selected-edge scores; retained transition-specific Cell controls",
        "redraw": "rendered from native phase arrays and joined score tables"}, indent=2) + "\n")


if __name__ == "__main__":
    main()
