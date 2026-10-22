#!/usr/bin/env python3
"""Display actual root-dependent DTI pseudotime and paired cell metrics.

Replots at full assembly width with equal-square maps, shared DPT scale,
and the same circle-to-diamond paired-score encoding and font sizes as C.
Renders saved paired refits and metric tables.
"""

import argparse
import hashlib
import json
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from plot_style import BLUE, GREY, style_axis, export


VARIANTS = ("original_root_refit", "hsc1_root_refit")
PANEL, GAP, MARGIN, FOOTER = 2.16, 0.16, 0.07, 0.48
ROW_WIDTH = 6.82
METRIC_LEFT = MARGIN + 2 * (PANEL + GAP) + 0.24
FONT = 6.7
METRICS = (("cbdir", "CBDir"), ("batc", "BATC"))


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def draw_dpt(ax, data, variant, limits):
    xy = np.asarray(data.obsm["X_umap"], dtype=float)
    values = data.obs["dpt_pseudotime"].to_numpy(dtype=float)
    order = np.argsort(values, kind="stable")
    ax.scatter(xy[order, 0], xy[order, 1], c=values[order], cmap="plasma",
               norm=Normalize(0, 1), s=1.7, alpha=0.96, linewidths=0,
               rasterized=True, zorder=1)
    labels = data.obs["clusters"].astype(str).to_numpy()
    offsets = {"HSC_1": (-11, 9), "Mono_1": (8, 10), "Ery_2": (9, -3)}
    for cluster, offset in offsets.items():
        center = np.median(xy[labels == cluster], axis=0)
        dx, dy = offset
        ax.annotate(cluster, center, xytext=offset, textcoords="offset points",
                    fontsize=FONT, color="#30343B", va="center",
                    ha="right" if dx < 0 else "left", zorder=4,
                    bbox=dict(facecolor="white", edgecolor="none", alpha=0.80, pad=0.25))
    root_index = int(data.uns["iroot"])
    marker, color = ("o", GREY) if variant == VARIANTS[0] else ("D", BLUE)
    ax.scatter(*xy[root_index], marker=marker, s=46, facecolor=color,
               edgecolor="white", linewidth=1.05, zorder=5)
    ax.set_xlim(*limits[0])
    ax.set_ylim(*limits[1])
    ax.set_aspect("equal", adjustable="box")
    ax.set_axis_off()


def draw_paired_metrics(ax, scores, summary):
    # Preserve the C panel's physical point separation on this narrower axis.
    # Offset the three seed pairs from the mean markers so all remain visible.
    jitter = np.array([-0.15, -0.09, 0.12])
    plotted = []
    for center, (metric, label) in zip((0, 1.15), METRICS):
        paired = scores.pivot(index="seed", columns="variant", values=metric).sort_index()
        paired = paired[list(VARIANTS)]
        if paired.shape != (3, 2) or not np.isfinite(paired.to_numpy()).all():
            raise ValueError(f"Expected exactly three complete seed pairs: {metric}")
        for index, seed in enumerate(paired.index):
            original, hsc = paired.loc[seed]
            ax.plot([center - 0.14 + jitter[index], center + 0.14 + jitter[index]],
                    [original, hsc], color="#A5ACB1", alpha=0.55, lw=0.5, zorder=1)
            plotted.append({"seed": int(seed), "metric": metric,
                            "original_root_refit": float(original), "hsc1_root_refit": float(hsc),
                            "paired_delta": float(hsc - original)})
        for index, (variant, color, marker) in enumerate(zip(VARIANTS, (GREY, BLUE), ("o", "D"))):
            x = center + (-0.14 if index == 0 else 0.14)
            values = paired[variant].to_numpy()
            ax.scatter(x + jitter, values, s=7, color=color, alpha=0.42,
                       edgecolor="none", marker=marker, zorder=2)
            mean = float(values.mean())
            expected = summary.set_index("metric").loc[metric, "original_mean" if index == 0 else "hsc1_mean"]
            np.testing.assert_allclose(mean, expected, atol=1e-12, rtol=0)
            ax.scatter(x, mean, s=40, color=color, marker=marker,
                       edgecolor="white", linewidth=0.7, zorder=4)
            ax.annotate(f"{mean:.2f}", (x, mean),
                        xytext=(-3 if index == 0 else 3, -8 if index == 0 else 9),
                        textcoords="offset points", ha="right" if index == 0 else "left",
                        va="center", color=color, fontsize=7)
    ax.axhline(0, color="#737B81", lw=0.7, zorder=0)
    ax.set(xlim=(-0.50, 1.65), ylim=(-0.08, 1.03), ylabel="Score")
    ax.set_xticks([0, 1.15], [item[1] for item in METRICS])
    ax.set_yticks([0, 0.5, 1])
    ax.grid(axis="y", color="#EEF0F2", lw=0.5)
    style_axis(ax)
    return plotted


def draw_guides(fig, width, height):
    cax = fig.add_axes([(MARGIN + 0.46) / width, 0.27 / height,
                        (2 * PANEL + GAP - 0.92) / width, 0.055 / height])
    bar = fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0, 1), cmap="plasma"),
                      cax=cax, orientation="horizontal", ticks=[0, 0.5, 1])
    bar.outline.set_visible(False)
    bar.ax.tick_params(labelsize=FONT, length=2, width=0.5, pad=2)
    bar.set_label("DPT", size=FONT, labelpad=2)
    handles = [Line2D([], [], ls="none", marker=marker, color=color,
                       markersize=4.5, label=label)
               for marker, color, label in (("o", GREY, "Original root"), ("D", BLUE, "HSC_1 root"))]
    fig.legend(handles=handles, frameon=False, loc="lower left",
               bbox_to_anchor=((METRIC_LEFT - 0.02) / width, 0.03 / height),
               fontsize=FONT, handletextpad=0.45, labelspacing=0.4, borderaxespad=0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20261001,
                        help="Paired-refit seed displayed in the maps; use the first recorded seed")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    provenance_path = args.data_root / "manifest.json"
    if not provenance_path.exists():
        provenance_path = args.data_root / "provenance.json"
    provenance = json.loads(provenance_path.read_text())
    if provenance["status"] != "complete" or args.seed != provenance["seeds"][0]:
        raise ValueError("Use completed experiment and the prespecified first paired seed")
    records, data = {}, {}
    for variant in VARIANTS:
        stem = f"{variant}_seed{args.seed}"
        record = provenance["fits"][stem]
        path = args.data_root / "native_outputs" / f"{stem}.h5ad"
        if sha256(path) != record["sha256"]:
            raise ValueError(f"Changed paired native result: {stem}")
        data[variant] = ad.read_h5ad(path)
        records[variant] = {"source": str(path), "sha256": record["sha256"], "root": record["root"]}
    left, right = (data[variant] for variant in VARIANTS)
    assert left.obs_names.equals(right.obs_names) and left.var_names.equals(right.var_names)
    for key in ("X_umap", "X_tsne", "X_diffmap"):
        np.testing.assert_array_equal(left.obsm[key], right.obsm[key])
    for key in ("palantir_pseudotime", "clusters"):
        np.testing.assert_array_equal(left.obs[key], right.obs[key])
    np.testing.assert_array_equal(left.obsm["X_umap"], left.obsm["X_tsne"])
    xy = np.asarray(left.obsm["X_umap"], dtype=float)
    if not np.isfinite(xy).all():
        raise ValueError("Invalid coordinates")
    for variant in VARIANTS:
        pseudotime = data[variant].obs["dpt_pseudotime"].to_numpy(dtype=float)
        if not np.isfinite(pseudotime).all() or pseudotime.min() < 0 or pseudotime.max() > 1:
            raise ValueError(f"Invalid native global DPT: {variant}")
        if pseudotime[int(data[variant].uns["iroot"])] != 0:
            raise ValueError(f"Native selected root has nonzero DPT: {variant}")
    center = (xy.min(axis=0) + xy.max(axis=0)) / 2
    half = float(np.max(xy.max(axis=0) - xy.min(axis=0)) * 0.55)
    limits = [(float(value - half), float(value + half)) for value in center]
    scores = pd.read_csv(args.data_root / "cell_scores.csv", float_precision="round_trip")
    scores = scores[scores.variant.isin(VARIANTS)].copy()
    if scores.duplicated(["variant", "seed"]).any() or len(scores) != 6:
        raise ValueError("Duplicate, missing, or extra fitted seed records")
    scores["seed"] = scores.seed.astype(int)
    if sorted(scores.seed.unique()) != sorted(provenance["seeds"]):
        raise ValueError("Seed set differs from paired experiment")
    summary = pd.read_csv(args.data_root / "paired_summary.csv", float_precision="round_trip")
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": FONT,
                         "axes.labelsize": 7.8, "xtick.labelsize": 6.8, "ytick.labelsize": 6.8,
                         "legend.fontsize": FONT, "pdf.fonttype": 42, "ps.fonttype": 42,
                         "svg.fonttype": "none", "axes.labelpad": 3})
    width, height = ROW_WIDTH, PANEL + FOOTER + MARGIN
    fig = plt.figure(figsize=(width, height), facecolor="white")
    for col, variant in enumerate(VARIANTS):
        ax = fig.add_axes([(MARGIN + col * (PANEL + GAP)) / width,
                           FOOTER / height, PANEL / width, PANEL / height])
        draw_dpt(ax, data[variant], variant, limits)
    ax = fig.add_axes([METRIC_LEFT / width, FOOTER / height,
                       (width - MARGIN - METRIC_LEFT) / width, PANEL / height])
    plotted = draw_paired_metrics(ax, scores, summary)
    draw_guides(fig, width, height)
    fig.canvas.draw()
    # As A/B: fixed physical canvas; no tight-bbox rescaling of the composite.
    for suffix in ("pdf", "svg", "png"):
        fig.savefig(args.output_dir / f"bonemarrow_root_control_row.{suffix}", dpi=300, facecolor="white")
    fig.savefig(args.output_dir / "bonemarrow_reading_size.png", dpi=150, facecolor="white")
    plt.close(fig)
    for variant in VARIANTS:
        single = plt.figure(figsize=(PANEL + 2 * MARGIN, PANEL + 2 * MARGIN), facecolor="white")
        ax = single.add_axes([MARGIN / (PANEL + 2 * MARGIN), MARGIN / (PANEL + 2 * MARGIN),
                              PANEL / (PANEL + 2 * MARGIN), PANEL / (PANEL + 2 * MARGIN)])
        draw_dpt(ax, data[variant], variant, limits)
        for suffix in ("pdf", "svg", "png"):
            single.savefig(args.output_dir / f"bonemarrow_{variant}.{suffix}", dpi=300, facecolor="white")
        plt.close(single)
    metric_fig, ax = plt.subplots(figsize=(2.15, 2.27))
    metric_fig.subplots_adjust(left=0.22, right=0.98, bottom=0.24, top=0.93)
    draw_paired_metrics(ax, scores, summary)
    export(metric_fig, args.output_dir / "bonemarrow_paired_cell_metrics")
    pd.DataFrame(plotted).to_csv(args.output_dir / "plotted_seed_pairs.csv", index=False)
    pd.DataFrame({"cell": left.obs_names, "tsne_1": xy[:, 0], "tsne_2": xy[:, 1],
                  "cluster": left.obs["clusters"].to_numpy(),
                  "palantir_pseudotime": left.obs["palantir_pseudotime"].to_numpy()}).to_csv(
                      args.output_dir / "plotted_cell_coordinates.csv", index=False)
    pd.DataFrame({"cell": left.obs_names,
                  **{variant: data[variant].obs["dpt_pseudotime"].to_numpy(dtype=float)
                     for variant in VARIANTS}}).to_csv(
                         args.output_dir / "plotted_dpt_values.csv", index=False)
    manifest = {
        "generator_sha256": sha256(Path(__file__)), "native_sources": records,
        "paired_fit_generator_sha256": provenance["script_sha256"], "displayed_seed": args.seed,
        "metric_seeds": provenance["seeds"], "n_cells": left.n_obs, "n_genes": left.n_vars,
        "plotting_style_source_sha256": sha256(Path(__file__).with_name("plot_style.py")),
        "map_style_source_sha256": sha256(Path(__file__).with_name("plot_tivelo_full_lineage_shared_dpt.py")),
        "source_hashes": {path.name: sha256(path) for path in
                          (provenance_path, args.data_root / "cell_scores.csv", args.data_root / "paired_summary.csv")},
        "coordinates": "archived t-SNE stored in X_umap; no new embedding, rotation or per-panel rescaling",
        "background": "Each condition's native obs.dpt_pseudotime; full-data DTI ordering, not path-specific orientation DPT",
        "dpt_column": "dpt_pseudotime", "shared_color_range": [0, 1],
        "dpt_renormalization": False, "streamlines": False,
        "dpt_value_hashes": {variant: hashlib.sha256(
            np.asarray(data[variant].obs["dpt_pseudotime"], dtype=float).tobytes()).hexdigest()
            for variant in VARIANTS},
        "pairing": "one-to-one optimizer seeds; not per-gene or biological replicates",
        "aggregation": "arithmetic mean of the three paired seeds; all raw seed points retained, no significance test",
        "displayed_cluster_labels": ["HSC_1", "Mono_1", "Ery_2"],
        "shared_limits": limits, "map_panel_inches": PANEL,
        "row_inches": [width, height], "gap_inches": GAP, "map_font_points": FONT,
        "metrics": [key for key, _ in METRICS], "baked_in_titles_captions_panel_letters": False,
    }
    (args.output_dir / "source_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2), flush=True)


if __name__ == "__main__":
    main()
