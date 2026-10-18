#!/usr/bin/env python3
"""Render gene comparison panels and summaries from registered inputs and CSVs."""

from __future__ import annotations

import argparse
import csv
import hashlib
from importlib.metadata import version
import json
import os
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import BoundaryNorm, ListedColormap, to_hex
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import sparse

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
from evaluation.gene.config import REAL_CONFIGS, SIM_CONFIGS, SIM_PATHS, VELOCITY_KEYS
from evaluation.gene.plotting import PHASE_PLOT_VERSION, plot_phase_portrait
from evaluation.cell.plots import BENCHMARK_COLORS

CASES = ["retina", "DentateGyrus_lamanno", "gastrulation", "capo_pre", "bonemarrow",
         "mousebrain", "pancreas", "branching", "murk", "all"]
# Join existing result-table identifiers to configuration IDs, not copied names.
CONFIG_IDS = {
    "scVelo_dyn": "C03", "MultiVelo": "C04", "UniTVelo": "C05",
    "TIvelo": "C08", "veloVI": "C12", "cellDancer": "C15",
    "DeepVelo": "C16", "DeepVelo_scv": "C18", "Graphvelo_scv": "C19",
}
METHODS = pd.read_csv(ROOT / "methods/configurations.tsv", sep="\t", keep_default_na=False)
ORDER = sorted(CONFIG_IDS, key=lambda method: METHODS.config_id.tolist().index(CONFIG_IDS[method]))
NAMES = {method: METHODS.set_index("config_id").loc[identifier, "manuscript_name"]
         for method, identifier in CONFIG_IDS.items()}
SUMMARY_NAMES = {**NAMES, "MultiVelo": "MultiVelo"}
COLORS = {SUMMARY_NAMES[method]: BENCHMARK_COLORS[NAMES[method]] for method in ORDER}
DATASETS = pd.read_csv(ROOT / "data/datasets.tsv", sep="\t", keep_default_na=False).set_index("dataset_id")
DATASET_IDS = {"branching": "sim_only_branching", "murk": "sim_only_murk", "all": "sim_hybrid"}
DATASET_NAMES = {case: DATASETS.loc[DATASET_IDS.get(case, case.lower()), "paper_name"] for case in CASES}
LAYOUTS = json.loads((HERE / "plotting/layouts.json").read_text())
CBDIR_CSV = HERE / "results/gene_cbdir_reference30_all.csv"
COSSIM_CSV = HERE / "results/simulation/gene_cbdir_cossim_reference.csv"


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_scores(cbdir_path=CBDIR_CSV, cossim_path=COSSIM_CSV):
    """Keep reference-30 CBDir and join frozen CosSim by the exact gene key."""
    scores = pd.read_csv(cbdir_path)
    required = {"dataset", "method", "gene", "new_cbdir", "source_h5ad", "velocity_s_key", "velocity_u_key"}
    if not required.issubset(scores):
        raise ValueError(f"CBDir table lacks columns: {sorted(required - set(scores))}")
    keys = ["dataset", "method", "gene"]
    if scores.duplicated(keys).any():
        raise ValueError("Duplicate dataset/method/gene keys in CBDir table")
    if set(scores.dataset) != set(CASES) or set(scores.method) != set(ORDER):
        raise ValueError("CBDir table does not match the registered gene cases/methods")
    sim = pd.read_csv(cossim_path).rename(columns={
        "Dataset": "dataset", "Method": "method", "Gene": "gene", "CosSim": "cossim",
    })
    if not set(keys + ["cossim"]).issubset(sim):
        raise ValueError("CosSim table lacks the gene keys or CosSim column")
    if sim.duplicated(keys).any():
        raise ValueError("Duplicate dataset/method/gene keys in CosSim table")
    expected = set(map(tuple, scores.loc[scores.dataset.isin(SIM_CONFIGS), keys].to_numpy()))
    if set(map(tuple, sim[keys].to_numpy())) != expected:
        raise ValueError("CosSim and reference-30 tables have different simulation gene keys")
    # The compact table's CBDir column is a different frozen result; do not
    # substitute it for the reference-30 new_cbdir used by these panels.
    scores = scores.merge(sim[keys + ["cossim"]], on=keys, how="left", validate="one_to_one")
    scores["display_method"] = scores.method.map(NAMES)
    scores["display_dataset"] = scores.dataset.map(DATASET_NAMES)
    return scores


def source_inputs(case, scores):
    """Check result locators against the shared evaluation input contract."""
    paths = REAL_CONFIGS[case]["paths"] if case in REAL_CONFIGS else SIM_PATHS[case]
    rows = []
    for method in ORDER:
        subset = scores[(scores.dataset == case) & (scores.method == method)]
        velocity_s, velocity_u = VELOCITY_KEYS[method]
        expected = {"source_h5ad": "artifact://velocity/" + paths[method],
                    "velocity_s_key": velocity_s, "velocity_u_key": velocity_u}
        for column, value in expected.items():
            recorded = subset[column].dropna().unique()
            if len(recorded) and set(recorded) != {value}:
                raise ValueError(f"{case}/{method}: {column} disagrees with evaluation config")
        rows.append({"method": method, "path": paths[method],
                     "velocity_s": velocity_s, "velocity_u": velocity_u})
    return rows


def read_phase_inputs(case, scores, artifact_root, genes):
    """Read selected coordinates into plotting copies, never rewriting H5ADs."""
    import anndata as ad

    reference = None
    if case in REAL_CONFIGS:
        reference = ad.read_h5ad(artifact_root / REAL_CONFIGS[case]["label_reference"], backed="r")
        try:
            labels = reference.obs["clusters"].copy()
            colors = reference.uns.get("clusters_colors")
        finally:
            reference.file.close()
    adatas, inputs = [], []
    for item in source_inputs(case, scores):
        path = artifact_root / item["path"]
        original = ad.read_h5ad(path)
        present = [gene for gene in genes if gene in original.var_names]
        obs = original.obs.copy()
        label_alignment = None
        if case in REAL_CONFIGS:
            if "Annotation" in obs:
                obs["clusters"] = obs["Annotation"].copy()
            aligned = labels.reindex(original.obs_names)
            if aligned.notna().all():
                obs["clusters"] = aligned.copy()
                label_alignment = "cell_id"
            elif len(labels) == original.n_obs:
                # Retained historical plotting fallback, recorded in provenance.
                obs["clusters"] = labels.copy()
                label_alignment = "equal_length_historical_fallback"
            else:
                raise ValueError(f"{case}/{item['method']}: reference labels cannot be aligned")
        elif "bin" not in obs and "clusters" in obs:
            obs["bin"] = obs["clusters"].copy()
        plot = ad.AnnData(X=np.zeros((original.n_obs, len(present)), dtype=np.float32),
                          obs=obs, var=pd.DataFrame(index=present))
        for target, source in [("Ms", "Ms"), ("Mu", "Mu"),
                               ("velocity", item["velocity_s"]), ("velocity_u", item["velocity_u"])]:
            if source not in original.layers:
                raise KeyError(f"{case}/{item['method']}: missing layer {source}")
            values = original[:, present].layers[source]
            values = values.toarray() if sparse.issparse(values) else np.asarray(values).copy()
            if case in REAL_CONFIGS and target.startswith("velocity"):
                # Match the old real-panel loader's display-only cleaning.
                values = np.clip(np.nan_to_num(values.astype(np.float64), nan=0.0,
                                              posinf=0.0, neginf=0.0), -1e10, 1e10)
            plot.layers[target] = values
        if case in REAL_CONFIGS and colors is not None:
            plot.uns["clusters_colors"] = list(colors)
        elif case in SIM_CONFIGS:
            # Retain the approved Figure 4 Early--Late progression palette.
            categories = list(plot.obs["bin"].cat.categories)
            stages = [int(str(label).rsplit("_", 1)[1]) for label in categories]
            if any(stage < 0 or stage > 9 for stage in stages):
                raise ValueError("Simulation progression colors require bins 0--9")
            palette = [to_hex(plt.get_cmap("plasma")(value)) for value in np.linspace(0.06, 0.94, 10)]
            plot.uns["bin_colors"] = [palette[stage] for stage in stages]
        inputs.append({"method": item["method"], "source": "artifact://velocity/" + item["path"],
                       "sha256": sha256(path), "velocity_s_key": item["velocity_s"],
                       "velocity_u_key": item["velocity_u"], "label_alignment": label_alignment})
        adatas.append(plot)
        del original
    return adatas, inputs


def plot_phase_panels(case, adatas, scores, genes, output_dir, dpi=600):
    import scvelo as scv

    scv.settings.set_figure_params("scvelo", dpi=dpi, frameon=False)
    plt.rcParams["figure.facecolor"] = "white"
    simulation = case in SIM_CONFIGS
    table = scores[scores.dataset == case].set_index(["method", "gene"])
    outputs = []
    for gene in genes:
        fig, axes = plt.subplots(1, len(ORDER), figsize=(45, 7.5))
        for ax, method, adata in zip(axes, ORDER, adatas):
            row = table.loc[(method, gene)]
            finite = np.isfinite(row.new_cbdir)
            metric = (f"$\\mathrm{{CBDir}}_{{\\mathrm{{gene}}}}$: {row.new_cbdir:.4f}"
                      if finite else "$\\mathrm{CBDir}_{\\mathrm{gene}}$: N/A")
            if simulation:
                metric += f"\nCosSim: {row.cossim:.4f}" if finite else "\nCosSim: N/A"
            if gene in adata.var_names:
                plot_phase_portrait(adata, gene, ax=ax, legend_loc=None,
                                    color_key="bin" if simulation else "clusters")
            ax.set_xlabel("")
            ax.set_ylabel("")
            ax.axis("off")
            ax.set_box_aspect(1)
            ax.set_title(NAMES[method], fontsize=20 if simulation else 28,
                         fontweight="bold", pad=45 if simulation else 55)
            ax.text(0.5, 1.02 if simulation else 1.04, metric, transform=ax.transAxes,
                    ha="center", va="bottom", fontsize=14 if simulation else 20)
        path = output_dir / f"{gene}_9methods_compare.png"
        plt.tight_layout()
        # Square phase axes settle their title positions during the first draw.
        # Validate the final text extents before exporting the fixed-size panel.
        for _ in range(3):
            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()
            highest = max(text.get_window_extent(renderer).y1
                          for ax in axes for text in [ax.title, *ax.texts])
            overflow = highest - fig.bbox.height + 12
            if overflow <= 0:
                break
            fig.subplots_adjust(top=fig.subplotpars.top - overflow / fig.bbox.height - 0.01)
        else:
            raise ValueError(f"Phase-panel titles exceed the canvas: {case}/{gene}")
        if case == "branching" and gene == "BoostALL_5":
            palette = [to_hex(plt.get_cmap("plasma")(value)) for value in np.linspace(0.06, 0.94, 10)]
            cmap = ListedColormap(palette)
            norm = BoundaryNorm(np.arange(-0.5, 10.5), cmap.N)
            cax = fig.add_axes([0.84, 0.055, 0.15, 0.018])
            colorbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax,
                                   orientation="horizontal", ticks=[-0.5, 9.5])
            colorbar.ax.set_xticklabels(["Early", "Late"])
            colorbar.ax.tick_params(labelsize=14, length=0, pad=3)
            colorbar.outline.set_visible(False)
            colorbar.ax.xaxis.set_label_position("top")
            colorbar.set_label("Cluster order", fontsize=14, labelpad=3)
        fig.savefig(path, dpi=dpi, facecolor="white", transparent=False)
        plt.close(fig)
        outputs.append(path)
    return outputs


def plot_real_summary(case, scores, output_path, dpi=600):
    frame = scores[(scores.dataset == case) & np.isfinite(scores.new_cbdir)].copy()
    frame["Method"] = frame.method.map(SUMMARY_NAMES)
    frame = frame.rename(columns={"gene": "Gene", "new_cbdir": "Overall_Score"})
    means = frame.groupby("Method").Overall_Score.mean()
    order = means.sort_values(ascending=False).index.tolist()
    heatmap = frame.pivot(index="Method", columns="Gene", values="Overall_Score").loc[order]
    layout = LAYOUTS[case]
    plt.rcParams.update(plt.rcParamsDefault)
    sns.set_theme(style="whitegrid", font_scale=1.1)
    fig, axes = plt.subplots(1, 2, figsize=(layout["fig_width"], layout["fig_height"]),
                            gridspec_kw={"width_ratios": [layout["boxplot_ratio"], layout["heatmap_ratio"]]})
    fig.patch.set_facecolor("white")
    sns.boxplot(x="Method", y="Overall_Score", data=frame, order=order, palette=COLORS, ax=axes[0])
    for index, method in enumerate(order):
        axes[0].scatter(index, means[method], color="#FFD700", s=48, marker="D",
                        edgecolor="black", linewidth=0.6, zorder=10)
    axes[0].set_title("")
    axes[0].set_xlabel("")
    axes[0].set_xticklabels(order, rotation=0, ha="center")
    axes[0].tick_params(axis="x", labelsize=layout["boxplot_label_size"], pad=10)
    sns.heatmap(heatmap, mask=heatmap.isna(), annot=True, cmap="RdBu_r", center=0,
                fmt=".2f", ax=axes[1], linewidths=0, linecolor="none",
                annot_kws={"size": layout["heatmap_annot_size"]}, cbar_kws={"label": ""})
    axes[1].set_facecolor("white")
    axes[1].grid(False)
    axes[1].legend(handles=[mpatches.Patch(facecolor="white", edgecolor="gray", label="Not Calculated")],
                   loc="lower right", bbox_to_anchor=(1.0, -0.33), frameon=False)
    for key in ("title", "ylabel", "xlabel"):
        getattr(axes[1], "set_" + key)("")
    axes[1].tick_params(axis="x", rotation=0, labelsize=layout["heatmap_label_size"], pad=8)
    plt.subplots_adjust(bottom=layout["bottom"], wspace=layout["wspace"])
    plt.tight_layout(rect=[0, 0.05, 1, 1])
    # The established export replaces the default ylabel after layout.
    axes[0].set_ylabel(r"$\mathrm{CBDir}_{\mathrm{gene}}$")
    fig.savefig(output_path, dpi=dpi, facecolor="white", transparent=False)
    plt.close(fig)


def plot_simulation_summary(case, scores, output_path, dpi=600):
    frame = scores[(scores.dataset == case) & np.isfinite(scores.new_cbdir)].copy()
    frame = frame.sort_values(["method", "gene"], key=lambda column: column.map(
        {name: i for i, name in enumerate(ORDER if column.name == "method" else SIM_CONFIGS[case]["genes"])}
    ))
    frame["Method"] = frame.method.map(SUMMARY_NAMES)
    frame = frame.rename(columns={"new_cbdir": "CBDir", "cossim": "CosSim"})
    sns.set(style="whitegrid", font_scale=1.2)
    fig, axes = plt.subplots(1, 2, figsize=(27, 7.5))
    for ax, metric in zip(axes, ["CBDir", "CosSim"]):
        means = frame.groupby("Method")[metric].mean()
        order = means.sort_values(ascending=False).index.tolist()
        sns.boxplot(x="Method", y=metric, data=frame, order=order, palette=COLORS,
                    showfliers=False, ax=ax)
        sns.stripplot(x="Method", y=metric, data=frame, order=order, color="black",
                      size=2, alpha=0.25, ax=ax)
        for index, method in enumerate(order):
            ax.scatter(index, means[method], color="#FFD700", s=48, marker="D",
                       edgecolor="black", linewidth=0.6, zorder=10)
        ax.set_title(r"$\mathrm{CBDir}_{\mathrm{gene}}$" if metric == "CBDir" else "CosSim")
        ax.set_xlabel("")
        ax.set_ylim(-1.0, 1.0)
        ax.tick_params(axis="x", rotation=0)
        if case == "murk" and metric == "CBDir":
            for tick in ax.get_xticklabels():
                if "_scv" in tick.get_text():
                    tick.set_fontsize(12)
    plt.tight_layout()
    fig.savefig(output_path, dpi=dpi, facecolor="white", transparent=False)
    plt.close(fig)


def write_manifest(output_root, cases, outputs):
    """Retain the existing paper-panel identities for generated assets."""
    rows = []
    mapping = {
        "Figure4_B": [("murk", gene) for gene in SIM_CONFIGS["murk"]["genes"]],
        "Figure4_C": [("branching", gene) for gene in SIM_CONFIGS["branching"]["genes"]],
        "Figure4_D": [(case, "summary") for case in SIM_CONFIGS],
        "Figure5_B": [("gastrulation", gene) for gene in REAL_CONFIGS["gastrulation"]["genes"]],
        "Figure5_C": [("gastrulation", "summary")],
        "Figure5_F": [("DentateGyrus_lamanno", gene) for gene in REAL_CONFIGS["DentateGyrus_lamanno"]["genes"]],
        "Figure5_G": [("retina", gene) for gene in REAL_CONFIGS["retina"]["genes"]],
        "Figure5_H": [(case, "summary") for case in ("retina", "DentateGyrus_lamanno")],
        "SuppFig1": [("DentateGyrus_lamanno", gene) for gene in REAL_CONFIGS["DentateGyrus_lamanno"]["genes"] if gene != "Ncald"],
        "SuppFig2": [("retina", gene) for gene in REAL_CONFIGS["retina"]["genes"] if gene != "Snap25"],
        "SuppFig3": [("gastrulation", gene) for gene in REAL_CONFIGS["gastrulation"]["genes"][1:11]],
        "SuppFig4": [("gastrulation", gene) for gene in REAL_CONFIGS["gastrulation"]["genes"][11:]],
        "SuppFig5": [("capo_pre", gene) for gene in REAL_CONFIGS["capo_pre"]["genes"]],
        "SuppFig6": [("bonemarrow", gene) for gene in REAL_CONFIGS["bonemarrow"]["genes"]],
        "SuppFig7": [("mousebrain", gene) for gene in REAL_CONFIGS["mousebrain"]["genes"]],
        "SuppFig8": [("pancreas", gene) for gene in REAL_CONFIGS["pancreas"]["genes"]],
    }
    for panel, assets in mapping.items():
        for case, asset in assets:
            if case not in cases:
                continue
            folder = output_root / ("real" if case in REAL_CONFIGS else "simulation") / case
            path = folder / ("summary_cbdir.png" if case in REAL_CONFIGS else "summary_cbdir_cossim.png") if asset == "summary" else folder / f"{asset}_9methods_compare.png"
            if path not in outputs:
                continue
            rows.append(dict(panel=panel, dataset=case, asset=asset,
                             path=path.relative_to(output_root).as_posix(),
                             score_source=CBDIR_CSV.relative_to(ROOT).as_posix(), sha256=sha256(path),
                             source_kind="supporting_gene_summary" if panel == "Figure4_D" else "score_annotated_panel"))
    fields = ["panel", "dataset", "asset", "path", "score_source", "sha256", "source_kind"]
    with (output_root / "assembly_manifest.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-root", type=Path, help="H5AD root, or VELOCITY_ARTIFACT_ROOT; unnecessary for summaries")
    parser.add_argument("--output-root", type=Path, default=ROOT / "outputs/gene-panels")
    parser.add_argument("--cases", nargs="+", choices=CASES, default=CASES)
    parser.add_argument("--genes", nargs="+", help="Selected registered genes for a focused phase-panel run")
    parser.add_argument("--summaries-only", action="store_true")
    parser.add_argument("--dpi", type=int, default=600)
    args = parser.parse_args()
    if args.dpi <= 0:
        parser.error("--dpi must be positive")
    if args.summaries_only and args.genes:
        parser.error("--genes requires phase panels; it does not restrict summary statistics")
    artifact = args.artifact_root or os.environ.get("VELOCITY_ARTIFACT_ROOT")
    if not args.summaries_only and not artifact:
        parser.error("phase panels require --artifact-root or VELOCITY_ARTIFACT_ROOT")
    output = args.output_root.expanduser().resolve()
    canonical = (HERE / "figures/panels").resolve()
    if output == canonical or canonical in output.parents:
        parser.error("render to a separate output directory; canonical paper assets require independent review")
    if output.exists() and any(output.iterdir()):
        parser.error("--output-root must be empty; use a new directory")
    genes_by_case = {case: (REAL_CONFIGS if case in REAL_CONFIGS else SIM_CONFIGS)[case]["genes"] for case in args.cases}
    if args.genes:
        for case, genes in genes_by_case.items():
            if not set(args.genes).issubset(genes):
                parser.error(f"--genes includes unregistered genes for {case}")
            genes_by_case[case] = args.genes
    scores = load_scores()
    output.mkdir(parents=True, exist_ok=True)
    outputs, inputs = [], {}
    for case in args.cases:
        folder = output / ("real" if case in REAL_CONFIGS else "simulation") / case
        folder.mkdir(parents=True, exist_ok=True)
        if not args.summaries_only:
            adatas, inputs[case] = read_phase_inputs(case, scores, Path(artifact).expanduser().resolve(), genes_by_case[case])
            outputs.extend(plot_phase_panels(case, adatas, scores, genes_by_case[case], folder, args.dpi))
            del adatas
        path = folder / ("summary_cbdir.png" if case in REAL_CONFIGS else "summary_cbdir_cossim.png")
        if case in REAL_CONFIGS:
            plot_real_summary(case, scores, path, args.dpi)
        else:
            plot_simulation_summary(case, scores, path, args.dpi)
        outputs.append(path)
        score_columns = ["dataset", "method", "gene", "new_cbdir", "cossim", "display_method", "display_dataset"]
        scores.loc[scores.dataset == case, score_columns].to_csv(folder / "scores_used.csv", index=False)
    write_manifest(output, args.cases, set(outputs))
    provenance = {
        "record_scope": "generation", "renderer": "evaluation/gene/plot_gene_assets.py",
        "renderer_sha256": sha256(__file__), "phase_plot_renderer": "evaluation/gene/plotting/phase_portrait.py",
        "phase_plot_renderer_sha256": sha256(HERE / "plotting/phase_portrait.py"),
        "phase_plot_version": PHASE_PLOT_VERSION, "cases": args.cases,
        "software_versions": {package: version(package) for package in
                              (["numpy", "pandas", "matplotlib", "seaborn", "scipy", "scikit-learn"]
                               + ([] if args.summaries_only else ["anndata", "scvelo"]))},
        "summaries_only": args.summaries_only, "dpi": args.dpi,
        "phase_genes": {} if args.summaries_only else genes_by_case,
        "score_sources": {p.relative_to(ROOT).as_posix(): sha256(p) for p in [CBDIR_CSV, COSSIM_CSV]},
        "inputs": inputs, "display_precision": {"gene_score": 4, "heatmap": 2},
        "real_phase_velocity_cleaning": "nan/inf to zero; clip to [-1e10, 1e10]; plotting copy only",
        "neighbors_recomputed": False, "moments_recomputed": False,
        "stripplot_jitter": "inherits NumPy display RNG; no effect on metric values or ordering",
        "outputs": {p.relative_to(output).as_posix(): sha256(p) for p in outputs},
    }
    for path in [HERE / "config.py", HERE / "plotting/layouts.json", ROOT / "methods/configurations.tsv", ROOT / "data/datasets.tsv", ROOT / "evaluation/cell/plots.py"]:
        provenance.setdefault("configuration_sources", {})[path.relative_to(ROOT).as_posix()] = sha256(path)
    (output / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(f"Rendered {len(outputs)} assets from registered CSVs and inputs into {output}")


if __name__ == "__main__":
    main()
