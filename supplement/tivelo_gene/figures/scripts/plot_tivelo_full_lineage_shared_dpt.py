#!/usr/bin/env python3
"""Export title-free TIvelo UMAP assets with complete lineage coverage.

All known developmental bins are colored. Main-path and secondary-path DPT
are shown separately, on identical original UMAP coordinates. No captions,
dataset headings or panel letters are baked into either assets or previews.
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
from matplotlib.patches import Polygon
import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull


BLUE, ORANGE, INK, GREY = "#0072B2", "#D55E00", "#30343B", "#D8DBDF"
SCENARIOS = (
    ("branching", "branching", "indicators/gene/simulation/branching_analyze/dataset/s1_branch.h5ad"),
    ("all", "hybrid", "indicators/gene/simulation/sim_all_analyze/data/sim_all.h5ad"),
)
KINDS = ("states", "known", "main_dpt", "branch_dpt")


def draw_umap(ax, item, kind):
    xy, labels = item["xy"], item["labels"]
    ax.scatter(xy[:, 0], xy[:, 1], s=1.5, color=GREY, alpha=0.90,
               linewidths=0, rasterized=True, zorder=1)
    if kind == "states":
        for node, color in ((item["early"], BLUE), (item["final"], ORANGE)):
            selected = xy[labels == node]
            ax.scatter(selected[:, 0], selected[:, 1], s=3.0, color=color,
                       linewidths=0, rasterized=True, zorder=3)
        selected = xy[labels == item["initial"]]
        center = selected.mean(axis=0)
        outline = selected[ConvexHull(selected).vertices]
        ax.add_patch(Polygon(center + 1.18 * (outline - center), closed=True,
                             fill=False, edgecolor=INK, linewidth=0.70,
                             joinstyle="round", zorder=4))
        offsets = (
            {item["early"]: (7, -1), item["initial"]: (6, 4), item["final"]: (-5, 9)}
            if item["scenario"] == "branching"
            else {item["early"]: (-9, 8), item["final"]: (0, 9)}
        )
        for node in dict.fromkeys((item["early"], item["initial"], item["final"])):
            color = BLUE if node == item["early"] else ORANGE if node == item["final"] else INK
            dx, dy = offsets[node]
            ax.annotate(node.replace("Branch", "").replace("_", ""),
                        xy[labels == node].mean(axis=0), xytext=(dx, dy),
                        textcoords="offset points", fontsize=6.7, color=color,
                        ha="right" if dx < 0 else "left", va="center", clip_on=False)
    else:
        if kind == "known":
            selected_xy, values = xy, item["bins"]
            norm = Normalize(0, 9)
        else:
            frame = item["frames"][kind]
            selected_xy = xy[item["positions"][kind]]
            values, norm = frame.dpt.to_numpy(dtype=float), Normalize(0, 1)
        order = np.argsort(values, kind="stable")
        ax.scatter(selected_xy[order, 0], selected_xy[order, 1], c=values[order],
                   cmap="plasma", norm=norm, s=1.7, alpha=0.96, linewidths=0,
                   rasterized=True, zorder=2)
    center, half = item["center"], item["half_range"]
    ax.set_xlim(center[0] - half, center[0] + half)
    ax.set_ylim(center[1] - half, center[1] + half)
    ax.set_aspect("equal", adjustable="box")
    ax.set_axis_off()


def render(items, kinds, output_stem, guides=True):
    panel, gap, margin = 1.58, 0.12, 0.07
    footer = 0.62 if guides else margin
    width = len(kinds) * panel + (len(kinds) - 1) * gap + 2 * margin
    height = len(items) * panel + (len(items) - 1) * 0.18 + footer + margin
    fig = plt.figure(figsize=(width, height), facecolor="white")
    for row, item in enumerate(items):
        y = footer + (len(items) - row - 1) * (panel + 0.18)
        for col, kind in enumerate(kinds):
            x = margin + col * (panel + gap)
            ax = fig.add_axes([x / width, y / height, panel / width, panel / height])
            draw_umap(ax, item, kind)

    if guides:
        for col, kind in enumerate(kinds):
            x = margin + col * (panel + gap)
            if kind == "states":
                handles = [
                    Line2D([], [], linestyle="none", marker="o", markersize=3.8,
                           markeredgewidth=0, markerfacecolor=BLUE, label="Early reference"),
                    Line2D([], [], linestyle="none", marker="o", markersize=4.7,
                           markeredgewidth=0.7, markeredgecolor=INK,
                           markerfacecolor="none", label="Initial origin"),
                    Line2D([], [], linestyle="none", marker="o", markersize=3.8,
                           markeredgewidth=0, markerfacecolor=ORANGE, label="Final root"),
                ]
                fig.legend(handles=handles, loc="lower left", frameon=False,
                           bbox_to_anchor=((x + 0.06) / width, 0.04 / height),
                           fontsize=6.7, handlelength=0.9, handletextpad=0.5,
                           labelspacing=0.4, borderaxespad=0)
            else:
                if kind == "branch_dpt" and "main_dpt" in kinds:
                    continue
                if kind == "main_dpt" and "branch_dpt" in kinds:
                    colorbar_left = x + 0.45
                    colorbar_width = 2 * panel + gap - 0.90
                else:
                    colorbar_left, colorbar_width = x + 0.18, 1.22
                cax = fig.add_axes([colorbar_left / width, 0.34 / height,
                                    colorbar_width / width, 0.065 / height])
                norm = Normalize(0, 9) if kind == "known" else Normalize(0, 1)
                ticks = [0, 3, 6, 9] if kind == "known" else [0, 0.5, 1]
                cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap="plasma"),
                                   cax=cax, orientation="horizontal", ticks=ticks)
                cbar.outline.set_visible(False)
                cbar.ax.tick_params(labelsize=6.7, length=2, width=0.5, pad=2)
                cbar.set_label("Developmental bin" if kind == "known" else "DPT",
                               size=6.7, labelpad=2)
    fig.canvas.draw()
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(output_stem.with_suffix(f".{suffix}"), dpi=300, facecolor="white")
    if len(kinds) == 4:
        fig.savefig(output_stem.parent / "ab_reading_size.png", dpi=120, facecolor="white")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-root", type=Path, required=True)
    parser.add_argument("--audit-root", type=Path, required=True)
    parser.add_argument("--branch-data", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    single_dir = args.output_dir / "individual"
    single_dir.mkdir(exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 6.7,
                         "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"})
    stages = pd.read_csv(args.audit_root / "stage_audit/stage_summary.csv").set_index(
        "scenario", verify_integrity=True)
    main_path = args.audit_root / "ordering_trigger_audit_20260929_threads1/cell_orders.csv"
    branch_path = args.branch_data / "branch_cell_orders.csv"
    tables = {"main_dpt": pd.read_csv(main_path, dtype={"cell": str}),
              "branch_dpt": pd.read_csv(branch_path, dtype={"cell": str})}
    items, source_rows = [], []
    for scenario, name, relative_path in SCENARIOS:
        adata = ad.read_h5ad(args.artifact_root / relative_path, backed="r")
        xy = np.asarray(adata.obsm["X_umap"], dtype=float).copy()
        labels = adata.obs["bin"].astype(str).to_numpy()
        cells = adata.obs_names.astype(str)
        adata.file.close()
        if not cells.is_unique or xy.shape != (len(cells), 2) or not np.isfinite(xy).all():
            raise ValueError(f"Invalid source UMAP or cell IDs: {scenario}")
        bins = np.array([int(label.rsplit("_", 1)[1]) for label in labels])
        frames, positions = {}, {}
        for kind, table in tables.items():
            path = "process_main" if kind == "main_dpt" else "branch_1"
            frame = table[(table.scenario == scenario) & (table.path == path)].copy()
            if frame.empty or frame.cell.duplicated().any():
                raise ValueError(f"Missing or duplicate source path cells: {scenario}/{kind}")
            index = cells.get_indexer(frame.cell)
            if np.any(index < 0) or not np.array_equal(bins[index], frame.bin.to_numpy()):
                raise ValueError(f"Source cell/bin mismatch: {scenario}/{kind}")
            source_branches = np.array([label.rsplit("_", 1)[0] for label in labels[index]])
            if not np.array_equal(source_branches, frame.branch.to_numpy()):
                raise ValueError(f"Source branch mismatch: {scenario}/{kind}")
            values = frame.dpt.to_numpy(dtype=float)
            if not np.isfinite(values).all() or values.min() < 0 or values.max() > 1:
                raise ValueError(f"Invalid DPT: {scenario}/{kind}")
            frames[kind], positions[kind] = frame, index
        if set(frames["main_dpt"].cell) | set(frames["branch_dpt"].cell) != set(cells):
            raise ValueError(f"Incomplete combined path coverage: {scenario}")
        stage = stages.loc[scenario]
        minimum, maximum = xy.min(axis=0), xy.max(axis=0)
        item = {"scenario": scenario, "name": name, "xy": xy, "labels": labels, "bins": bins,
                "early": "BranchA_0", "initial": str(stage.process_start_node),
                "final": str(stage.correct_start_node), "frames": frames, "positions": positions,
                "center": (minimum + maximum) / 2,
                "half_range": float(np.max(maximum - minimum) * 0.55)}
        if item["final"] != stage.dti_start_node:
            raise ValueError(f"Final root differs from DTI root: {scenario}")
        items.append(item)
        source_rows.append({"scenario": scenario, "source_h5ad": relative_path,
                            "known_order_colored_cells": len(cells),
                            "main_dpt_colored_cells": len(frames["main_dpt"]),
                            "branch_dpt_colored_cells": len(frames["branch_dpt"]),
                            "shared_path_cells": len(set(frames["main_dpt"].cell) & set(frames["branch_dpt"].cell)),
                            "umap_float64_sha256": hashlib.sha256(xy.tobytes()).hexdigest(),
                            "original_main_dpt_sha256": hashlib.sha256(frames["main_dpt"].dpt.to_numpy(dtype=float).tobytes()).hexdigest()})

    render(items, ("states",), args.output_dir / "state_selection_material")
    render(items, ("known", "main_dpt", "branch_dpt"), args.output_dir / "order_comparison_material")
    render(items, KINDS, args.output_dir / "ab_material_preview")
    for item in items:
        for kind in KINDS:
            render([item], (kind,), single_dir / f"{item['name']}_{kind}", guides=False)

    manifest = {"purpose": "Title-free rendering assets with full lineage coverage",
                "generator_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                "main_order_csv_sha256": hashlib.sha256(main_path.read_bytes()).hexdigest(),
                "branch_order_csv_sha256": hashlib.sha256(branch_path.read_bytes()).hexdigest(),
                "datasets": source_rows, "preview_rows": [item["name"] for item in items],
                "preview_columns": list(KINDS), "independent_path_dpt": True,
                "known_order": "All cells, discrete developmental bin 0-9; not exact chronological time",
                "grey_cells_in_dpt": "Cells outside the independently analyzed path",
                "shared_path_dpt_colorbar": True,
                "umap_transform": "None; identical square limits within a dataset",
                "baked_in_titles_captions_panel_letters": False}
    (args.output_dir / "source_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(source_rows, indent=2), flush=True)


if __name__ == "__main__":
    main()
