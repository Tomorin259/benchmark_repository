"""Render Supplementary Figure S15 materials and Supplementary Table S9.

Reads committed compact diagnostic inputs and exports panels and table rows.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import sys

import anndata as ad
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.collections import PathCollection
from matplotlib.patches import Circle

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))
from evaluation.gene.plotting.phase_portrait import plot_phase_portrait
from evaluation.cell.plots import BENCHMARK_COLORS
from summarize_controls import manuscript_table


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def style(ax):
    ax.spines[["top", "right"]].set_visible(False)
    for side in ["left", "bottom"]:
        ax.spines[side].set_color("#9EA3A7")
        ax.spines[side].set_linewidth(.65)
    ax.tick_params(length=2.6, width=.65, pad=3)


def export(fig, folder, stem, dpi):
    for extension in ["png", "pdf", "svg"]:
        target = folder / f"{stem}.{extension}"
        fig.savefig(target, dpi=dpi, facecolor="white")
        if extension == "svg":
            target.write_text("\n".join(line.rstrip() for line in target.read_text().splitlines()) + "\n")
    plt.close(fig)


def plot_paired_controls(root, output, dpi):
    # Use the accepted S11 PCA layout: neutral labels, compact filled/open
    # pairs, faint method-colored links, black means and one shared key.
    refits = pd.read_csv(root / "results/multivelo_refits.csv")
    scv = pd.read_csv(root / "results/scvelo_controls.csv")
    names = pd.read_csv(REPO / "methods/configurations.tsv", sep="\t").set_index("operational_id").manuscript_name
    width_inches = 7.2
    fig = plt.figure(figsize=(width_inches, width_inches * 420 / 1060))
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set(xlim=(0, 1060), ylim=(420, 0))
    ax.axis("off")
    pixel_to_pt = width_inches * 72 / 1060

    def line(x, y, u, v, color="#333333", width=1, opacity=1):
        ax.plot((x, u), (y, v), color=color, lw=width * pixel_to_pt,
                alpha=opacity, solid_capstyle="butt")

    def text(x, y, label, anchor="center", size=16, rotation=0):
        ax.text(x, y, label, ha=anchor, va="center" if rotation else "baseline",
                color="#262626", fontsize=size * pixel_to_pt,
                fontfamily="DejaVu Sans", rotation=rotation)

    def dot(x, y, color, hollow=False, radius=3.7):
        ax.add_patch(Circle((x, y), radius, facecolor="white" if hollow else color,
                            edgecolor=color if hollow else "#333333", lw=pixel_to_pt))

    plotted = []
    for j, (dataset, title) in enumerate([("branching", "normal + branching"), ("all", "hybrid: branching")]):
        left, top, bottom = 70 + j * 530, 40, 325
        right = left + 435
        y = lambda value: bottom - float(value) / 1.03 * (bottom - top)
        for tick in [0, .2, .4, .6, .8, 1]:
            line(left, y(tick), right, y(tick), "#E7E7E7", .7)
            text(left - 10, y(tick) + 5, f"{tick:.1f}", "right", 14)
        line(left, top, left, bottom, width=1.4)
        line(left, bottom, right, bottom, width=1.4)
        text((left + right) / 2, 23, title, size=18)
        for i, (method, name, color, frame, field, original, changed) in enumerate([
            ("scVelo_dyn", names["scvelo_dynamic"], BENCHMARK_COLORS[names["scvelo_dynamic"]], scv, "variant", "saved", "force_induction_fit_tau_smoothed"),
            ("MultiVelo", names["multivelo"], BENCHMARK_COLORS[names["multivelo"]], refits, "condition", "refit_default", "on_allow_repression"),
        ]):
            center = left + 110 + i * 215
            frame = frame[(frame.dataset == dataset) & ~frame.gene.str.startswith("OneSidedBoost")]
            original_values = frame[frame[field] == original].set_index("gene").sort_index()
            changed_values = frame[frame[field] == changed].set_index("gene").sort_index()
            assert len(original_values) == 20 and original_values.index.equals(changed_values.index)
            for k, gene in enumerate(original_values.index):
                a, b = center - 24 + k - 9.5, center + 24 + k - 9.5
                first = float(original_values.loc[gene, "cossim"])
                second = float(changed_values.loc[gene, "cossim"])
                assert 0 <= first <= 1 and 0 <= second <= 1
                line(a, y(first), b, y(second), color, .8, .35)
                dot(a, y(first), color)
                dot(b, y(second), color, True)
            for x, data, condition, anchor, offset in [
                (center - 24, original_values, original, "right", -16),
                (center + 24, changed_values, changed, "left", 16),
            ]:
                mean = float(data.cossim.mean())
                line(x - 12, y(mean), x + 12, y(mean), "#262626", 2)
                text(x + offset, y(mean) + 4, f"{mean:.4f}", anchor, 13)
                for gene, row in data.iterrows():
                    plotted.append(dict(dataset=dataset, method=method, condition=condition, gene=gene, cbdir=float(row.cbdir), cossim=float(row.cossim)))
            text(center, bottom + 26, name, size=14)
    text(18, 185, "CosSim", rotation=90)
    dot(845, 390, "#666666")
    text(858, 395, "Default", "left", 14)
    dot(960, 390, "#666666", True)
    text(1035, 395, "Modified", "right", 14)
    assert all(item.get_color() == "#262626" for item in ax.texts)
    export(fig, output, "paired_controls", dpi)
    pd.DataFrame(plotted).to_csv(output / "scores_used.csv", index=False)


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--module-root", type=Path, default=HERE)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--dpi", type=int, default=600)
    parser.add_argument("--paired-only", action="store_true", help="Render the paired-score row without touching phase portraits or tables")
    args = parser.parse_args()
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise ValueError("Use an empty output directory")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    root = args.module_root
    inputs = [root / "inputs/phase_case.npz", root / "results/multivelo_refits.csv", root / "results/scvelo_controls.csv", root / "results/condition_summary.csv", root / "results/candidate_gate_effects.csv", root / "CAPTION.md"]
    before = {str(p.relative_to(root)): sha(p) for p in inputs}
    font = "Arial" if "Arial" in {f.name for f in font_manager.fontManager.ttflist} else "DejaVu Sans"
    plt.rcParams.update({"font.family": font, "font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.linewidth": .65, "pdf.fonttype": 42, "svg.fonttype": "none"})
    if args.paired_only:
        used = {name: value for name, value in before.items() if name in ["results/multivelo_refits.csv", "results/scvelo_controls.csv"]}
        plot_paired_controls(root, args.output_dir, args.dpi)
        assert all(sha(root / name) == value for name, value in before.items())
        record = dict(inputs=used, input_hashes_unchanged=True, script_sha256=sha(__file__), font="DejaVu Sans", dpi=args.dpi, panel_letters=False, primary_results_modified=False, text_color="#262626", style_reference=dict(path="supplement/veloae/plot_pca_summary.py", sha256=sha(REPO / "supplement/veloae/plot_pca_summary.py")), method_registry_sha256=sha(REPO / "methods/configurations.tsv"), condition_key={"Default": "saved native output or matching default refit", "Modified": {"scVelo_dyn": "force_induction_fit_tau_smoothed", "MultiVelo": "on_allow_repression"}}, outputs={p.name: sha(p) for p in args.output_dir.iterdir() if p.is_file()})
        (args.output_dir / "manifest.json").write_text(json.dumps(record, indent=2) + "\n")
        print("Rendered paired-score materials:", args.output_dir)
        return
    with np.load(root / "inputs/phase_case.npz", allow_pickle=False) as stored:
        arrays = {k: stored[k] for k in stored.files}
    obj = ad.AnnData(X=np.zeros_like(arrays["Ms"]))
    obj.obs_names = arrays["cell_ids"].astype(str)
    obj.var_names = ["case_gene"]
    for key in ["Ms", "Mu", "velocity_s_true", "velocity_u_true", "default_velo_s", "default_velo_u", "gate_open_velo_s", "gate_open_velo_u"]:
        obj.layers[key] = arrays[key]
    obj.obs["known_progression"] = pd.Categorical(arrays["bins"], categories=list(range(10)))
    obj.uns["known_progression_colors"] = [plt.get_cmap("plasma")(i / 9) for i in range(10)]
    refits = pd.read_csv(root / "results/multivelo_refits.csv")
    case_scores = refits[(refits.dataset == "branching") & (refits.gene == "BoostALL_5")].set_index("condition")
    diagnostics = []
    phase_specs = [
        ("phase_truth", "Ground truth\n", "velocity_s_true", "velocity_u_true", None),
        ("phase_multivelo_default", "MultiVelo (RNA_only)\nDefault", "default_velo_s", "default_velo_u", "refit_default"),
        ("phase_multivelo_gate_open", "MultiVelo (RNA_only)\nGate open", "gate_open_velo_s", "gate_open_velo_u", "on_allow_repression"),
    ]
    for stem, title, vs, vu, condition in phase_specs:
        fig = plt.figure(figsize=(2.35, 2.65))
        ax = fig.add_axes([.18, .20, .78, .64])
        _, diagnostic = plot_phase_portrait(obj, "case_gene", ax=ax, layer_v_s=vs, layer_v_u=vu, color_key="known_progression", legend_loc=None, s=2.7, alpha=.50, return_diagnostics=True)
        for collection in list(ax.collections):
            if isinstance(collection, PathCollection) and len(collection.get_facecolors()) == 0:
                collection.remove()
        ax.set_title(title, fontsize=8, pad=8, linespacing=1.45)
        ax.set_xlabel("Spliced (s)")
        ax.set_ylabel("Unspliced (u)" if condition is None else "")
        if condition is not None:
            ax.tick_params(labelleft=False)
            ax.text(.98, .035, f"CosSim {case_scores.loc[condition, 'cossim']:.4f}", transform=ax.transAxes, ha="right", va="bottom", fontsize=7)
        style(ax)
        diagnostic["gene"] = "BoostALL_5"
        diagnostics.append(diagnostic)
        export(fig, args.output_dir, stem, args.dpi)
    assert diagnostics[0]["sampled_cell_ids"] == diagnostics[1]["sampled_cell_ids"] == diagnostics[2]["sampled_cell_ids"]

    fig = plt.figure(figsize=(1.175, .42))
    bar = fig.add_axes([.05, .53, .90, .17])
    cb = fig.colorbar(ScalarMappable(norm=Normalize(0, 9), cmap="plasma"), cax=bar, orientation="horizontal", ticks=[0, 9])
    cb.ax.set_xticklabels(["Early", "Late"])
    cb.ax.get_xticklabels()[0].set_ha("left")
    cb.ax.get_xticklabels()[1].set_ha("right")
    cb.ax.tick_params(length=0, pad=2, labelsize=7)
    cb.outline.set_visible(False)
    export(fig, args.output_dir, "progression_key", args.dpi)

    plot_paired_controls(root, args.output_dir, args.dpi)

    (args.output_dir / "tableS9.tex").write_text(manuscript_table(root))
    assert all(sha(root / name) == value for name, value in before.items())
    manifest = dict(inputs=before, input_hashes_unchanged=True, script_sha256=sha(__file__), font=font, dpi=args.dpi, phase_diagnostics=diagnostics, panel_letters=False, primary_results_modified=False, figure="S15", table="S9", phase_score_scope="Mean of cell-wise CosSim over all 2000 cells of BoostALL_5; branch-specific scores use 1000 cells", style_reference=dict(path="supplement/veloae/plot_pca_summary.py", sha256=sha(REPO / "supplement/veloae/plot_pca_summary.py")), method_registry_sha256=sha(REPO / "methods/configurations.tsv"), condition_key={"Default": "saved native output or matching default refit", "Modified": {"scVelo_dyn": "force_induction_fit_tau_smoothed", "MultiVelo": "on_allow_repression"}}, outputs={p.name: sha(p) for p in args.output_dir.iterdir() if p.is_file()})
    (args.output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print("Rendered S15 materials, shared key, and Table S9:", args.output_dir)


if __name__ == "__main__":
    main()
