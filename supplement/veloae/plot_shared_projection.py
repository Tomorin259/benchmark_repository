"""Render the VeloAE PCA/velocity projection panels and their 5-column plate.

The NPZ files are trusted exports from the VeloAE analysis.  In particular,
``X2`` and ``V2`` are displayed as supplied; this script does not refit PCA.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator, ScalarFormatter


ALIASES = {
    "DentateGyrus_lamanno": "dentategyrus_lamanno",
    "bonemarrow": "bonemarrow",
    "capofibro": "capo_pre",
    "dentategyrus": "dentategyrus",
    "forebrain": "forebrain",
    "gastrulation": "gastrulation",
    "hematopoiesis": "hematopoiesis",
    "mousebrain": "mousebrain",
    "neuron": "scnt_neuron",
    "organoid": "sceu_organoid",
    "pancreas": "pancreas",
    "retina": "retina",
    "rpe1": "sceu_rpe1",
    "rpe1_gen": "rpe1_gen",
}
ORDER = tuple(ALIASES)
COLUMNS, ROWS = 5, 3
PANEL_WIDTH_MM, PANEL_HEIGHT_MM, FOOTER_MM = 45, 43, 9


def _registry_names(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    by_id = {row["dataset_id"]: row["paper_name"] for row in rows}
    names = {}
    for stem, dataset_id in ALIASES.items():
        if dataset_id not in by_id:
            raise ValueError(f"dataset {dataset_id!r} missing from registry {path}")
        names[stem] = by_id[dataset_id]
    return names


def _find_npz(input_dir: Path, stem: str, dataset_id: str) -> Path:
    direct = input_dir / f"{stem}.npz"
    if not direct.exists():
        raise FileNotFoundError(f"required NPZ is missing: {direct}")
    return direct


def _colors(categories: np.ndarray) -> dict[str, tuple[float, float, float, float]]:
    values = [str(x) for x in categories.tolist()]
    if len(set(values)) != len(values):
        raise ValueError("label_categories contains duplicates")
    if len(values) <= 20:
        palette = plt.get_cmap("tab20")(np.arange(len(values)))
    else:
        palette = plt.get_cmap("hsv")(np.linspace(0.0, 1.0, len(values), endpoint=False))
    return dict(zip(values, palette))


def _arrow_bins(xy: np.ndarray, velocity: np.ndarray, xlim, ylim) -> tuple[np.ndarray, np.ndarray, int]:
    n = xy.shape[0]
    minimum = max(5, n // 450)
    x_edges = np.linspace(*xlim, 15)
    y_edges = np.linspace(*ylim, 11)
    bins = []
    for iy in range(10):
        for ix in range(14):
            mask = (
                (xy[:, 0] >= x_edges[ix])
                & (xy[:, 0] <= x_edges[ix + 1] if ix == 13 else xy[:, 0] < x_edges[ix + 1])
                & (xy[:, 1] >= y_edges[iy])
                & (xy[:, 1] <= y_edges[iy + 1] if iy == 9 else xy[:, 1] < y_edges[iy + 1])
            )
            if int(mask.sum()) >= minimum:
                bins.append((float(xy[mask, 0].mean()), float(xy[mask, 1].mean()),
                             float(velocity[mask, 0].mean()), float(velocity[mask, 1].mean()),
                             iy, ix, int(mask.sum())))
    bins.sort(key=lambda row: (row[1], row[0]))
    if len(bins) > 28:
        picks = np.floor(np.linspace(0, len(bins)-1, 28)).astype(int)
        bins = [bins[i] for i in picks]
    points = np.asarray([row[:2] for row in bins], dtype=float).reshape(-1, 2)
    vectors = np.asarray([row[2:4] for row in bins], dtype=float).reshape(-1, 2)
    return points, vectors, minimum


def _formatter(ax: plt.Axes) -> None:
    xfmt = ScalarFormatter(useOffset=False, useMathText=False)
    yfmt = ScalarFormatter(useOffset=False, useMathText=False)
    for fmt in (xfmt, yfmt):
        fmt.set_scientific(True)
        fmt.set_powerlimits((-3, 4))
    ax.xaxis.set_major_formatter(xfmt)
    ax.yaxis.set_major_formatter(yfmt)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
    ax.tick_params(labelsize=6.5, length=2, pad=1.5)
    ax.xaxis.get_offset_text().set_fontsize(6.5)
    ax.yaxis.get_offset_text().set_fontsize(6.5)


def _plot_panel(stem: str, paper_name: str, data: dict[str, np.ndarray], out: Path) -> dict:
    xy = np.asarray(data["X2"], dtype=float)
    velocity = np.asarray(data["V2"], dtype=float)
    labels = np.asarray(data["labels"], dtype=str)
    categories = np.asarray(data["label_categories"], dtype=str)
    if xy.ndim != 2 or xy.shape[1] != 2 or velocity.shape != xy.shape:
        raise ValueError(f"{stem}: X2/V2 must both have shape (n, 2)")
    if labels.shape != (xy.shape[0],):
        raise ValueError(f"{stem}: labels length does not match X2")
    if not np.isfinite(xy).all() or not np.isfinite(velocity).all():
        raise ValueError(f"{stem}: X2/V2 contain non-finite values")
    if len(categories) < 2 or any(not str(value) for value in categories):
        raise ValueError(f"{stem}: label_categories must contain at least two non-empty categories")
    if not set(labels).issubset(set(categories)):
        raise ValueError(f"{stem}: labels include categories absent from label_categories")
    if len(set(labels)) < 2:
        raise ValueError(f"{stem}: only one observed annotation category")
    colors = _colors(categories)
    xspan = max(float(np.ptp(xy[:, 0])), np.finfo(float).eps)
    yspan = max(float(np.ptp(xy[:, 1])), np.finfo(float).eps)
    xmin, xmax = float(xy[:, 0].min()), float(xy[:, 0].max())
    ymin, ymax = float(xy[:, 1].min()), float(xy[:, 1].max())
    # independent bounds: 4% horizontal and 8% vertical padding
    xlim = (xmin - 0.04 * xspan, xmax + 0.04 * xspan)
    ylim = (ymin - 0.08 * yspan, ymax + 0.08 * yspan)
    fig, ax = plt.subplots(figsize=(45 / 25.4, 43 / 25.4), dpi=300)
    for category in categories:
        mask = labels == category
        ax.scatter(xy[mask, 0], xy[mask, 1], s=2.0, c=[colors[str(category)]], alpha=0.42,
                   linewidths=0, rasterized=True)
    arrow_xy, arrow_v, minimum = _arrow_bins(xy, velocity, xlim, ylim)
    magnitude = np.linalg.norm(arrow_v, axis=1)
    positive = magnitude[np.isfinite(magnitude) & (magnitude > 0)]
    p80 = float(np.percentile(positive, 80)) if positive.size else 1.0
    scale = 0.58 * (xspan / 14.0) / max(p80, np.finfo(float).eps)
    shown = arrow_v * scale
    ax.set_xlim(*xlim); ax.set_ylim(*ylim)
    ax.set_xlabel("", fontsize=7, labelpad=1)
    ax.set_ylabel("", fontsize=7, labelpad=1)
    title = paper_name.replace('_', '_\n', 1) if len(paper_name) > 24 and '_' in paper_name else '\n'.join(textwrap.wrap(paper_name, width=24))
    ax.set_title(title, fontsize=7.5, pad=2)
    _formatter(ax)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.subplots_adjust(left=0.18, right=0.98, bottom=0.17, top=0.86)
    # Cap arrows in rendered axes pixels, preserving their display direction.
    transform = ax.transData
    origin = transform.transform(arrow_xy)
    endpoint = transform.transform(arrow_xy + shown)
    display_vector = endpoint - origin
    display_length = np.linalg.norm(display_vector, axis=1)
    max_pixels = 0.12 * ax.get_window_extent().width
    factors = np.minimum(1.0, max_pixels / np.maximum(display_length, np.finfo(float).eps))
    endpoint = origin + display_vector * factors[:, None]
    shown = transform.inverted().transform(endpoint) - arrow_xy
    ax.quiver(arrow_xy[:, 0], arrow_xy[:, 1], shown[:, 0], shown[:, 1], angles="xy",
              scale_units="xy", scale=1, color="#303030", width=0.0032,
              headwidth=3.8, headlength=4.5, headaxislength=4.0, zorder=4)
    png = out / f"panel_{ORDER.index(stem) + 1:02d}_{stem}.png"
    fig.savefig(png, dpi=300, facecolor="white")
    fig.savefig(out / f"panel_{ORDER.index(stem) + 1:02d}_{stem}.pdf", facecolor="white")
    fig.savefig(out / f"panel_{ORDER.index(stem) + 1:02d}_{stem}.svg", facecolor="white")
    plt.close(fig)
    return {
        "dataset_id": ALIASES[stem],
        "paper_name": paper_name,
        "n_cells": int(xy.shape[0]),
        "n_categories": int(len(categories)),
        "categories": categories.tolist(),
        "x2_range": [float(xy[:, 0].min()), float(xy[:, 0].max())],
        "y2_range": [float(xy[:, 1].min()), float(xy[:, 1].max())],
        "v2_range": [[float(velocity[:, i].min()), float(velocity[:, i].max())] for i in range(2)],
        "explained_variance_ratio": [float(v) for v in np.asarray(data["explained_variance_ratio"]).ravel()[:2]],
        "arrow_count": int(len(arrow_xy)),
        "minimum_bin_count": int(minimum),
        "arrow_scale": float(scale),
        "axis_limits": [list(xlim), list(ylim)],
        "grid_bins": [14, 10],
        "arrow_display_cap_axes_width": 0.12,
        "png_dpi": 300,
        "panel_size_mm": [45, 43],
        "npz_sha256": "",
    }


def _assemble(out: Path) -> None:
    paths = [out / f"panel_{i:02d}_{stem}.png" for i, stem in enumerate(ORDER, 1)]
    from PIL import Image
    panel_images = [Image.open(path).convert("RGB") for path in paths]
    width, height = panel_images[0].size
    # Shared key and one orientation schematic; no per-panel repetition.
    for name, size in (('shared_key', (45, 43)), ('orientation', (45, 9))):
        fig = plt.figure(figsize=(size[0]/25.4, size[1]/25.4), dpi=300)
        ax = fig.add_axes((0, 0, 1, 1))
        ax.set(xlim=(0, size[0]), ylim=(0, size[1])); ax.axis('off')
        props = dict(arrowstyle='-|>', lw=.7, color='#303030', mutation_scale=6,
                     shrinkA=0, shrinkB=0)
        if name == 'shared_key':
            ax.scatter([9], [27], s=10, color='#888888', linewidths=0)
            ax.text(14, 27, 'Cells', fontsize=7, va='center')
            ax.annotate('', xy=(12, 19), xytext=(6, 19), arrowprops=props)
            ax.text(14, 19, 'Projected\nvelocity', fontsize=7, va='center')
        else:
            ax.annotate('', xy=(10, 2.2), xytext=(4, 2.2), arrowprops=props)
            ax.annotate('', xy=(4, 7.8), xytext=(4, 2.2), arrowprops=props)
            ax.text(11, 2.2, 'PC1', fontsize=7, va='center')
            ax.text(5, 7.2, 'PC2', fontsize=7, va='center')
        for suffix in ('png', 'pdf', 'svg'):
            fig.savefig(out/f'{name}.{suffix}', dpi=300, facecolor='white')
        plt.close(fig)
    orientation = Image.open(out/'orientation.png').convert('RGB')
    key = Image.open(out/'shared_key.png').convert('RGB')
    plate = Image.new("RGB", (width * COLUMNS, height * ROWS + orientation.height), "white")
    for i, image in enumerate(panel_images):
        plate.paste(image, ((i % COLUMNS) * width, (i // COLUMNS) * height))
    plate.paste(key, ((COLUMNS-1)*width, (ROWS-1)*height))
    plate.paste(orientation, (0, ROWS*height))
    plate.save(out / "veloae_pca_velocity_5x3.png", dpi=(300, 300))
    for image in panel_images:
        image.close()
    fig = plt.figure(figsize=(COLUMNS*PANEL_WIDTH_MM/25.4,
                              (ROWS*PANEL_HEIGHT_MM+FOOTER_MM)/25.4), dpi=300)
    ax = fig.add_axes((0, 0, 1, 1)); ax.imshow(plate); ax.axis('off')
    fig.savefig(out / "veloae_pca_velocity_5x3.pdf", facecolor="white")
    fig.savefig(out / "veloae_pca_velocity_5x3.svg", facecolor="white")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    names = _registry_names(args.registry)
    manifest = {"registry_name": args.registry.name, "input_name": args.input_dir.name, "panels": []}
    for stem in ORDER:
        path = _find_npz(args.input_dir, stem, ALIASES[stem])
        with np.load(path, allow_pickle=True) as loaded:
            required = {"X2", "V2", "labels", "label_categories", "W", "explained_variance_ratio"}
            missing = required - set(loaded.files)
            if missing:
                raise KeyError(f"{path}: missing keys {sorted(missing)}")
            data = {key: loaded[key] for key in loaded.files}
        record = _plot_panel(stem, names[stem], data, args.output_dir)
        record["npz_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
        manifest["panels"].append(record)
    _assemble(args.output_dir)
    (args.output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
