#!/usr/bin/env python3
"""Render Cell comparison streams from Source Data arrays and native outputs."""
import argparse
import base64
import gc
import io
import hashlib
import json
import sys
from pathlib import Path

import anndata as ad
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont
from matplotlib import font_manager
import scvelo as scv

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
CASES = {
    "gastrulation": ("real", "A", "Figure_2", "gastrulation_20methods_3x7_BATC"),
    "DentateGyrus_lamanno": ("real", "B", "Figure_2", "DentateGyrus_lamanno_20methods_3x7_CBDir"),
    "murk": ("simulation", "A", "Supplementary_Figure_12", "A_murk_streams"),
    "branching": ("simulation", "B", "Supplementary_Figure_12", "B_branching_streams"),
    "all": ("simulation", "C", "Supplementary_Figure_12", "C_all_streams"),
}


def draw_stream(ax, arrays, velocity, categories, colors):
    coordinates = arrays["X_umap"]
    obj = ad.AnnData(np.zeros((len(coordinates), 1), dtype=np.float32))
    obj.obs["clusters"] = pd.Categorical(arrays["group"].astype(str), categories=categories, ordered=True)
    obj.uns["clusters_colors"] = colors
    obj.obsm["X_umap"] = coordinates
    scv.pl.velocity_embedding_stream(
        obj, basis="umap", X=coordinates, V=np.nan_to_num(velocity),
        color="clusters", palette=colors, title="", legend_loc="none",
        frameon=False, show=False, ax=ax, density=2,
    )
    # Masked stream-grid endpoints are omitted consistently in raster and vector exports.
    for collection in ax.collections:
        widths = np.asarray(collection.get_linewidths())
        if not np.isfinite(widths).all():
            collection.set_linewidths(np.where(np.isfinite(widths), widths, 0.0))
    for patch in ax.patches:
        if not np.isfinite(patch.get_linewidth()):
            patch.set_visible(False)
    spread = np.ptp(coordinates, axis=0)
    ax.set_xlim(coordinates[:, 0].min() - .05 * spread[0], coordinates[:, 0].max() + .05 * spread[0])
    ax.set_ylim(coordinates[:, 1].min() - .05 * spread[1], coordinates[:, 1].max() + .05 * spread[1])
    ax.set_axis_off()


# Accepted Figure 2 montage geometry, before uniform 1/3 export scaling.
REAL_LAYOUTS = {
    "gastrulation": {"width": 24750, "height": 11139, "title": 874, "panel": 2771},
    "DentateGyrus_lamanno": {"width": 24690, "height": 11112, "title": 872, "panel": 2764},
}
REAL_STREAM_CROP = (612, 180, 3492, 2076)
REAL_ROW_GAP = 102
REAL_EXPORT_SCALE = 3


def fit_font(draw, text, font_path, width, maximum):
    font = ImageFont.truetype(font_path, maximum)
    text_width = draw.textbbox((0, 0), text, font=font)[2]
    size = max(1, int(maximum * min(1, .9 * width / max(text_width, 1))))
    return ImageFont.truetype(font_path, size)


def render_real_montage(dataset, inputs, scores, arrays, velocity, output):
    """Use accepted common crops and tile geometry for the real-data examples."""
    from evaluation.cell.plots import METHOD_NAME_MAPPING
    layout = REAL_LAYOUTS[dataset]
    base_width, remainder = divmod(layout["width"], 7)
    widths = [base_width + int(i < remainder) for i in range(7)]
    positions = np.cumsum([0, *widths[:-1]]).tolist()
    title_height, panel_height = layout["title"], layout["panel"]
    cell_height = title_height + panel_height
    scale = base_width / round((REAL_STREAM_CROP[2] - REAL_STREAM_CROP[0]) * .84)
    font_bold = font_manager.findfont(font_manager.FontProperties(family="Liberation Sans", weight="bold"), fallback_to_default=False)
    font_regular = font_manager.findfont(font_manager.FontProperties(family="Liberation Sans"), fallback_to_default=False)
    probe = ImageDraw.Draw(Image.new("RGB", (1, 1)))
    method_names = [METHOD_NAME_MAPPING[r.operational_id] for r in inputs.itertuples(index=False)]
    title_font = fit_font(probe, max(method_names, key=len), font_bold, min(widths), round(330 * scale * 1.15))
    metric_font = fit_font(probe, "CBDir: -0.0000    ICCoh: 0.0000", font_regular, min(widths), min(round(230 * scale), round(title_font.size * .82)))
    # The native panel follows the same 6x4 source frame and shared crop as
    # the retained individual-method panels, rather than a montage subplot.
    scv.set_figure_params("scvelo")
    obj = ad.AnnData(np.zeros((len(arrays["cell_ids"]), 1), dtype=np.float32))
    palettes = json.loads((ROOT / "evaluation/cell/stream_preview_config.json").read_text())[dataset]
    obj.obs["clusters"] = pd.Categorical(arrays["group"].astype(str), categories=palettes["categories"], ordered=True)
    obj.uns["clusters_colors"] = palettes["palette"]
    obj.obsm["X_umap"] = arrays["X_umap"]
    ax = scv.pl.velocity_embedding_stream(obj, basis="umap", X=arrays["X_umap"], V=velocity,
        color="clusters", palette=palettes["palette"], title="", legend_loc="none",
        frameon=False, dpi=600, show=False)
    buffer = io.BytesIO()
    ax.figure.savefig(buffer, format="png", dpi=600)
    plt.close(ax.figure)
    buffer.seek(0)
    native = Image.open(buffer)
    native.load()
    if native.size != (3600, 2400):
        raise ValueError(f"Native stream source frame differs: {native.size}")
    canvas = Image.new("RGB", (layout["width"] // REAL_EXPORT_SCALE, layout["height"] // REAL_EXPORT_SCALE), "white")
    reference_path = ROOT / f"evaluation/cell/figures/references/{dataset}.png"
    with Image.open(reference_path) as reference:
        if reference.size != (widths[0], cell_height):
            raise ValueError("Reference panel must use the complete accepted first tile")
        canvas.paste(reference.resize((round(widths[0] / 3), round(cell_height / 3)), Image.Resampling.LANCZOS), (0, 0))
    source_hashes = {str(reference_path.relative_to(ROOT)): hashlib.sha256(reference_path.read_bytes()).hexdigest()}
    for index, row in enumerate(inputs.itertuples(index=False), 1):
        tile_row, column = divmod(index, 7)
        width = widths[column]
        tile = Image.new("RGB", (width, cell_height), "white")
        draw = ImageDraw.Draw(tile)
        score = scores.loc[dataset, row.operational_id]
        texts = [METHOD_NAME_MAPPING[row.operational_id], f"CBDir: {score.cbdir:.4f}    ICCoh: {score.iccoh:.4f}", f"BATC: {score.batc:.4f}"]
        fonts = [title_font, metric_font, metric_font]
        boxes = [draw.textbbox((0, 0), text, font=font) for text, font in zip(texts, fonts)]
        heights = [b[3] - b[1] for b in boxes]
        gap = round(18 * scale * .75)
        metric_gap = max(gap + 30, round(gap * 2.8))
        top = max(18, gap, title_height - sum(heights) - gap - metric_gap - max(18, gap))
        for i, (text, font, box) in enumerate(zip(texts, fonts, boxes)):
            draw.text(((width - box[2] + box[0]) // 2, top - box[1]), text, font=font, fill=(20, 20, 20) if i == 0 else (70, 76, 85))
            top += heights[i] + (gap if i == 0 else metric_gap)
        source = native if row.operational_id == "multivelo" else Image.open(ROOT / f"evaluation/cell/streams/{dataset}/{row.operational_id}.png")
        if source.size != (3600, 2400):
            raise ValueError(f"Unexpected retained source frame: {row.operational_id}/{source.size}")
        cropped = source.crop(REAL_STREAM_CROP).convert("RGBA")
        background = Image.new("RGBA", cropped.size, "white")
        background.alpha_composite(cropped)
        panel = background.resize((width, panel_height), Image.Resampling.LANCZOS).convert("RGB")
        tile.paste(panel, (0, title_height))
        canvas.paste(tile.resize((round(width / 3), round(cell_height / 3)), Image.Resampling.LANCZOS),
                     (round(positions[column] / 3), round(tile_row * (cell_height + REAL_ROW_GAP) / 3)))
        if row.operational_id != "multivelo":
            source_path = ROOT / f"evaluation/cell/streams/{dataset}/{row.operational_id}.png"
            source_hashes[str(source_path.relative_to(ROOT))] = hashlib.sha256(source_path.read_bytes()).hexdigest()
            source.close()
        del tile, cropped, background, panel
    canvas.save(output.with_suffix(".png"), dpi=(200, 200))
    canvas.save(output.with_suffix(".pdf"), resolution=200)
    encoded = base64.b64encode(output.with_suffix(".png").read_bytes()).decode()
    output.with_suffix(".svg").write_text(f'<svg xmlns="http://www.w3.org/2000/svg" width="{canvas.width}px" height="{canvas.height}px" viewBox="0 0 {canvas.width} {canvas.height}"><image width="{canvas.width}" height="{canvas.height}" href="data:image/png;base64,{encoded}"/></svg>\n')
    native.close()
    del canvas, obj
    gc.collect()
    return source_hashes


def native_projection(row, template, artifact_root):
    relative = row.result_path.removeprefix("velocity:")
    if not row.result_path.startswith("velocity:") or Path(relative).is_absolute() or ".." in Path(relative).parts:
        raise ValueError("Expected a safe velocity: output locator")
    with h5py.File(artifact_root / relative) as handle:
        obs = handle["obs"]
        names = obs[obs.attrs["_index"]].asstr()[:]
        lookup = {name: i for i, name in enumerate(names)}
        positions = [lookup[str(name)] for name in template["cell_ids"]]
        coordinates = handle["obsm"][row.source_embedding_key][:][positions, :2]
        velocity = handle["obsm"][row.velocity_embedding_key][:][positions, :2]
    if not np.allclose(coordinates, template["X_umap"], rtol=0, atol=1e-5):
        raise ValueError("Native projection and published display coordinates differ")
    return velocity


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-data-root", type=Path, required=True)
    parser.add_argument("--artifact-root", type=Path, required=True, help="The velocity artifact directory")
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--dataset", action="append", choices=sorted(CASES))
    args = parser.parse_args()
    if args.output_root.exists() and any(args.output_root.iterdir()):
        parser.error("Use an empty output directory")
    args.output_root.mkdir(parents=True, exist_ok=True)
    registry = pd.read_csv(ROOT / "methods/configurations.tsv", sep="\t").set_index("config_id")
    palettes = json.loads((ROOT / "evaluation/cell/stream_preview_config.json").read_text())
    scv.settings.verbosity = 0
    plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42, "svg.fonttype": "none"})
    rendered = []
    sources = {}
    retained_panels = {}
    for dataset in args.dataset or list(CASES):
        scope, letter, folder, stem = CASES[dataset]
        module = ROOT / "evaluation/cell" / ("simulation" if scope == "simulation" else "")
        inputs = pd.read_csv(module / "inputs.tsv", sep="\t")
        inputs = inputs[inputs.evaluation_dataset.eq(dataset)].sort_values("config_id")
        scores = pd.read_csv(module / "results/cell_metrics.csv").set_index(["dataset", "method"])
        if scope == "real":
            row = next(r for r in inputs.itertuples(index=False) if r.operational_id == "multivelo")
            source = args.source_data_root / "main_figures/Figure_2" / f"{letter}_{dataset}_C04_stream.npz"
            sources[str(source.relative_to(args.source_data_root))] = hashlib.sha256(source.read_bytes()).hexdigest()
            with np.load(source, allow_pickle=False) as loaded:
                arrays = {key: loaded[key] for key in loaded.files}
            velocity = native_projection(row, arrays, args.artifact_root)
            retained_panels.update(render_real_montage(dataset, inputs, scores, arrays, velocity, args.output_root / stem))
            rendered.extend(stem + "." + suffix for suffix in ["png", "pdf", "svg"])
            continue
        fig, axes = plt.subplots(1, 9, figsize=(18, 2.4))
        fig.subplots_adjust(left=.003, right=.997, bottom=.02, top=.78, wspace=.08)
        slots = list(axes)
        if len(inputs) != len(slots):
            raise ValueError(f"Incomplete comparison: {dataset}")
        for ax, row in zip(slots, inputs.itertuples(index=False)):
            subfolder = "main_figures" if scope == "real" else "supplementary_figures"
            source = args.source_data_root / subfolder / folder / f"{letter}_{dataset}_{row.config_id}_stream.npz"
            sources[str(source.relative_to(args.source_data_root))] = hashlib.sha256(source.read_bytes()).hexdigest()
            with np.load(source, allow_pickle=False) as loaded:
                arrays = {key: loaded[key] for key in loaded.files}
            if row.operational_id == "multivelo":
                velocity = native_projection(row, arrays, args.artifact_root)
            elif scope == "simulation":
                velocity = arrays["velocity_umap_saved"]
            else:
                velocity = arrays.get("velocity_umap_display", arrays.get("velocity_umap_native"))
            if scope == "simulation":
                categories = arrays["group_categories"].astype(str).tolist()
                colors = arrays["group_colors"].astype(str).tolist()
            else:
                categories, colors = palettes[dataset]["categories"], palettes[dataset]["palette"]
            draw_stream(ax, arrays, velocity, categories, colors)
            score = scores.loc[dataset, row.operational_id]
            title = registry.loc[row.config_id, "manuscript_name"]
            if scope == "real":
                title += f"\nCBDir: {score.cbdir:.4f}   ICCoh: {score.iccoh:.4f}\nBATC: {score.batc:.4f}"
            else:
                title += f"\nCBDir: {score.cbdir:.4f}\nBATC: {score.batc:.4f}"
            ax.set_title(title, fontsize=10 if scope == "real" else 7)
        for suffix in ["png", "pdf", "svg"]:
            path = args.output_root / (stem + "." + suffix)
            fig.savefig(path, dpi=300, facecolor="white")
            rendered.append(path.name)
        plt.close(fig)
    provenance = {
        "generator": "evaluation/cell/plot_cell_stream_panels.py",
        "generator_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "source_data_arrays": sources,
        "retained_real_stream_panels": retained_panels,
        "real_montage_geometry": REAL_LAYOUTS,
        "real_shared_crop_pixels": list(REAL_STREAM_CROP),
        "real_uniform_export_scale": 1 / REAL_EXPORT_SCALE,
        "native_output_selection": "evaluation/cell/native_outputs.tsv",
        "native_output_selection_sha256": hashlib.sha256((ROOT / "evaluation/cell/native_outputs.tsv").read_bytes()).hexdigest(),
        "score_tables": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                         for path in [ROOT / "evaluation/cell/results/cell_metrics.csv",
                                      ROOT / "evaluation/cell/simulation/results/cell_metrics.csv"]},
        "outputs": {name: hashlib.sha256((args.output_root / name).read_bytes()).hexdigest() for name in rendered},
    }
    (args.output_root / "generation.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print("Rendered", len(rendered), "Cell stream material files")


if __name__ == "__main__":
    main()
