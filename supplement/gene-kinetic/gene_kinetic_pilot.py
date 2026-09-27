#!/usr/bin/env python3
"""Render masked nine-method gene spliced--unspliced phase portraits."""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import scanpy as sc


REPO = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
ARTIFACT_ROOT = None
REAL_PATHS = {}
DIAPH3_SCORES = None
METHODS = [
    ("cellDancer", "celldancer", "celldancer"),
    ("DeepVelo", "deepvelo", "deepvelo"),
    ("DeepVelo_scv", "deepvelo2022", "deepvelo2022"),
    ("GraphVelo_scv", "graphvelo_scv", "graphvelo_scv"),
    ("MultiVelo", "multivelo", "multivelo"),
    ("scVelo_dynamic", "scvelo_dyn", "scvelo_dynamic"),
    ("TIvelo", "tivelo", "tivelo"),
    ("UniTVelo", "unitvelo_with_vu", "unitvelo"),
    ("veloVI", "veloVI", "veloVI"),
]
CASES = {
    "gastrulation_hbax": {
        "dataset": "Gastrulation", "token": "Gastrulation", "gene": "Hba-x",
        "clusters": ["Erythroid1", "Erythroid2", "Erythroid3"],
        "edge": ["Erythroid2", "Erythroid3"], "score_edges": [["Erythroid1", "Erythroid2"], ["Erythroid2", "Erythroid3"]],
    },
    "dentate_ncald": {
        "dataset": "DentateGyrus_lamanno", "token": "DentateGyrus_lamanno", "gene": "Ncald",
        "clusters": ["Nbl2", "ImmGranule1", "ImmGranule2", "Granule"],
        "edge": ["ImmGranule1", "ImmGranule2"], "score_edges": [["ImmGranule1", "ImmGranule2"], ["ImmGranule2", "Granule"]],
    },
    "capo_e2f7": {
        "dataset": "capo_pre", "token": "capo_pre", "gene": "E2F7",
        "clusters": ["G1", "S", "G2M"], "edge": ["G1", "S"], "score_edges": [["G1", "S"], ["S", "G2M"]],
    },
}

# Additional README-prescribed gene/edge cases. Existing three pilot cases are
# retained under their original names so their output paths remain stable.
CASES.update({
    "gastrulation_gypa": {"dataset":"Gastrulation", "token":"Gastrulation", "gene":"Gypa", "clusters":["Erythroid1","Erythroid2","Erythroid3"], "edge":["Erythroid2","Erythroid3"], "score_edges":[["Erythroid1","Erythroid2"],["Erythroid2","Erythroid3"]]},
    "gastrulation_alas2": {"dataset":"Gastrulation", "token":"Gastrulation", "gene":"Alas2", "clusters":["Erythroid1","Erythroid2","Erythroid3"], "edge":["Erythroid2","Erythroid3"], "score_edges":[["Erythroid1","Erythroid2"],["Erythroid2","Erythroid3"]]},
    "gastrulation_fam210b": {"dataset":"Gastrulation", "token":"Gastrulation", "gene":"Fam210b", "clusters":["Erythroid1","Erythroid2","Erythroid3"], "edge":["Erythroid2","Erythroid3"], "score_edges":[["Erythroid1","Erythroid2"],["Erythroid2","Erythroid3"]]},
    "gastrulation_klf3": {"dataset":"Gastrulation", "token":"Gastrulation", "gene":"Klf3", "clusters":["Erythroid1","Erythroid2","Erythroid3"], "edge":["Erythroid2","Erythroid3"], "score_edges":[["Erythroid1","Erythroid2"],["Erythroid2","Erythroid3"]]},
    "gastrulation_hemgn": {"dataset":"Gastrulation", "token":"Gastrulation", "gene":"Hemgn", "clusters":["Erythroid1","Erythroid2","Erythroid3"], "edge":["Erythroid2","Erythroid3"], "score_edges":[["Erythroid1","Erythroid2"],["Erythroid2","Erythroid3"]]},
    "gastrulation_fech": {"dataset":"Gastrulation", "token":"Gastrulation", "gene":"Fech", "clusters":["Erythroid1","Erythroid2","Erythroid3"], "edge":["Erythroid2","Erythroid3"], "score_edges":[["Erythroid1","Erythroid2"],["Erythroid2","Erythroid3"]]},
    "dentate_diaph3": {"dataset":"DentateGyrus_lamanno", "token":"DentateGyrus_lamanno", "gene":"Diaph3", "clusters":["nIPC","Nbl1","Nbl2"], "edge":["nIPC","Nbl1"], "score_edges":[["nIPC","Nbl1"],["Nbl1","Nbl2"]]},
    "dentate_dab1": {"dataset":"DentateGyrus_lamanno", "token":"DentateGyrus_lamanno", "gene":"Dab1", "clusters":["CA","CA1-Sub","CA2-3-4"], "edge":["CA","CA1-Sub"], "score_edges":[["CA","CA1-Sub"]]},
    "capo_fam111a": {"dataset":"capo_pre", "token":"capo_pre", "gene":"FAM111A", "clusters":["G1","S","G2M"], "edge":["G1","S"], "score_edges":[["G1","S"],["S","G2M"]]},
    "retina_celf4": {"dataset":"retina", "token":"retina", "gene":"Celf4", "clusters":["Neuroblast","RGC","AC/HC"], "edge":["Neuroblast","RGC"], "score_edges":[["Neuroblast","RGC"],["Neuroblast","AC/HC"]]},
    "mousebrain_flrt2": {"dataset":"mousebrain", "token":"mousebrain", "gene":"Flrt2", "clusters":["V-SVZ","Upper Layer","Deeper Layer"], "edge":["V-SVZ","Upper Layer"], "score_edges":[["V-SVZ","Upper Layer"]]},
    "pancreas_g6pc2": {"dataset":"Pancreas", "token":"Pancreas", "gene":"G6pc2", "clusters":["Pre-endocrine","Beta","Alpha"], "edge":["Pre-endocrine","Beta"], "score_edges":[["Pre-endocrine","Beta"],["Pre-endocrine","Alpha"]]},
})
EDGE_SCORES = None

def load_edge_scores():
    path = REPO / "evaluation/gene/results/gene_cbdir_edge_details_reference30_all.csv"
    return pd.read_csv(path), path

def edge_labels(dataset, method, gene, edges):
    rows = EDGE_SCORES[(EDGE_SCORES.dataset == dataset) & (EDGE_SCORES.method == method) & (EDGE_SCORES.gene == gene)]
    labels = []
    for source, target in edges:
        row = rows[(rows.source == source) & (rows.target == target)]
        if gene == "Diaph3" and source == "nIPC" and target == "Nbl1" and DIAPH3_SCORES is not None:
            row = DIAPH3_SCORES[(DIAPH3_SCORES.method == method) & (DIAPH3_SCORES.source == source) & (DIAPH3_SCORES.target == target) & DIAPH3_SCORES.status.eq("ok")]
            score_column = "edge_score" if gene == "Diaph3" and source == "nIPC" else "new_edge_score"
        else:
            score_column = "new_edge_score"
        labels.append(f"{source}→{target}: " + ("%+.3f" % float(row.iloc[0][score_column]) if not row.empty and np.isfinite(row.iloc[0][score_column]) else "N/A"))
    return labels


def method_paths(token):
    dataset_key = {"Gastrulation":"gastrulation", "DentateGyrus_lamanno":"DentateGyrus_lamanno", "Pancreas":"pancreas", "capo_pre":"capo_pre", "retina":"retina", "mousebrain":"mousebrain"}.get(token)
    if dataset_key in REAL_PATHS:
        return {display: "velocity/" + REAL_PATHS[dataset_key][registry_name] for display, registry_name in [("cellDancer","cellDancer"),("DeepVelo","DeepVelo"),("DeepVelo_scv","DeepVelo_scv"),("GraphVelo_scv","Graphvelo_scv"),("MultiVelo","MultiVelo"),("scVelo_dynamic","scVelo_dyn"),("TIvelo","TIvelo"),("UniTVelo","UniTVelo"),("veloVI","veloVI")]}
    common = {
        "cellDancer": f"velocity/results/celldancer/celldancer_{token}.h5ad",
        "DeepVelo": f"velocity/results/deepvelo/results_with_vu/deepvelo_{token}.h5ad",
        "DeepVelo_scv": f"velocity/results/deepvelo2022/deepvelo_{token}_recon.h5ad",
        "GraphVelo_scv": f"velocity/results/graphvelo_scv/graphvelo{token}.h5ad",
        "MultiVelo": f"velocity/results/multivelo/multivelo_{token}.h5ad",
        "scVelo_dynamic": f"velocity/results/scvelo_dynamic/scv_dynamic_{token}.h5ad",
        "TIvelo": f"velocity/results/TIvelo/tivelo_{token}.h5ad",
        "UniTVelo": f"velocity/results/unitvelo/results_with_vu/unitvelo_{token}_with_velocity_u_mode1.h5ad",
        "veloVI": f"velocity/results/veloVI/veloVI_{token}_all.h5ad",
    }
    # These are the established historical exceptions used by Fig. 3/4.
    if token == "Gastrulation":
        common["DeepVelo_scv"] = "velocity/results/deepvelo2022/deepvelo_Gastrulation_recon.h5ad"
        common["veloVI"] = "velocity/results/veloVI/veloVI_Gastrulation_all.h5ad"
    if token == "retina":
        common["DeepVelo_scv"] = "velocity/results/deepvelo2022/deepvelo_retina_full_recon.h5ad"
        common["veloVI"] = "velocity/results/veloVI/veloVI_retina_full_all.h5ad"
    if token == "DentateGyrus_lamanno":
        common["DeepVelo_scv"] = "velocity/results/deepvelo2022/deepvelo_DentateGyrus_lamanno_recon.h5ad"
        common["veloVI"] = "velocity/results/veloVI/veloVI_DentateGyrus_lamanno_all.h5ad"
    return common


VELOCITY_KEYS = {
    "cellDancer": ("velocity", "velocity_u"), "DeepVelo": ("velocity", "velocity_unspliced"),
    "DeepVelo_scv": ("velocity_reconstructed", "velocity_u"), "GraphVelo_scv": ("velocity_gv", "velocity_u"),
    "MultiVelo": ("velo_s", "velo_u"), "scVelo_dynamic": ("velocity", "velocity_u"),
    "TIvelo": ("velocity", "velocity_u"), "UniTVelo": ("velocity", "velocity_u"), "veloVI": ("velocity", "velocity_u"),
}

def align_layers(adata, ref_clusters):
    if "Annotation" in adata.obs:
        adata.obs["clusters"] = adata.obs["Annotation"]
    aligned = ref_clusters.reindex(adata.obs_names)
    if aligned.notna().all():
        adata.obs["clusters"] = aligned
    # Keep method-native matrices sparse until the selected-cluster phase
    # subset is made; phase_portrait densifies only the requested gene vector.
    return adata

def masked_gene_adata(adata, ref_clusters, display, gene, selected):
    source_s, source_u = VELOCITY_KEYS[display]
    if gene not in adata.var_names or source_s not in adata.layers or source_u not in adata.layers:
        return None, {"status": "unavailable"}
    mask = np.asarray(adata.obs["clusters"].astype(str).isin(selected))
    sub = adata[mask, [gene]].copy()
    diagnostics = {"status": "ok", "n_obs": int(sub.n_obs)}
    for out_key, source_key in [("velocity", source_s), ("velocity_u", source_u)]:
        x = sub.layers[source_key]
        x = x.toarray() if hasattr(x, "toarray") else np.asarray(x)
        x = np.asarray(x, dtype=np.float64)
        nonfinite = ~np.isfinite(x)
        diagnostics[f"{out_key}_nonfinite_raw"] = int(nonfinite.sum())
        diagnostics[f"{out_key}_replaced_zero"] = int(nonfinite.sum())
        x = np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)
        x = np.clip(x, -1e10, 1e10)
        sub.layers[out_key] = x
    for key in ["Ms", "Mu"]:
        x = sub.layers[key]
        sub.layers[key] = x.toarray() if hasattr(x, "toarray") else np.asarray(x)
    return sub, diagnostics


def make_panel(case_name, case, module):
    os.chdir(ARTIFACT_ROOT.parent)
    paths = method_paths(case["token"])
    ref_path = paths["DeepVelo"]
    ref = sc.read_h5ad(ref_path)
    ref_clusters = ref.obs["clusters"].copy()
    ref_colors = ref.uns.get("clusters_colors")
    adatas, source_paths = [], []
    preprocessing_diagnostics = {}
    for display, _, _ in METHODS:
        path = Path(paths[display])
        source_path = (ARTIFACT_ROOT.parent / path).resolve()
        source_paths.append(source_path)
        if not source_path.is_file():
            adatas.append(None)
            continue
        a = align_layers(sc.read_h5ad(source_path), ref_clusters)
        if ref_colors is not None:
            a.uns["clusters_colors"] = list(ref_colors)
        sub, prep_diag = masked_gene_adata(a, ref_clusters, display, case["gene"], case["clusters"])
        preprocessing_diagnostics[display] = prep_diag
        del a
        adatas.append(sub)
    valid = [a for a in adatas if a is not None]
    if not valid:
        raise RuntimeError(f"No method files found for {case_name}")
    selected = case["clusters"]
    out = OUT / "figures"; out.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 5, figsize=(14, 6.5), constrained_layout=True)
    phase_mod = module
    score_manifest = {}
    phase_diagnostics = {}
    dataset_score = {"Gastrulation":"gastrulation", "Pancreas":"pancreas"}.get(case["dataset"], case["dataset"])
    for ax, (display, _, _), a in zip(axes.flat[:9], METHODS, adatas):
        score_method = {"GraphVelo_scv":"Graphvelo_scv", "scVelo_dynamic":"scVelo_dyn"}.get(display, display)
        labels = edge_labels(dataset_score, score_method, case["gene"], case["score_edges"])
        score_manifest[display] = labels
        ax.set_title(display + "\n" + "\n".join(labels), fontsize=8, pad=2)
        if a is None or not all(k in a.layers for k in ["Ms", "Mu", "velocity", "velocity_u"]):
            ax.text(.5, .5, "gene unavailable", ha="center", va="center", transform=ax.transAxes); ax.axis("off"); continue
        sub = a
        _, diagnostics = phase_mod.plot_phase_portrait(sub, case["gene"], ax=ax, color_key="clusters", legend_loc=None, s=12, alpha=.60, arrow_grid=(25,25), arrow_density_percentile=15.0, arrow_length_mode="log1p", arrow_reference_quantile=.90, arrow_reference_length=.05, arrow_max_length=.10, return_diagnostics=True)
        phase_diagnostics[display] = {key: diagnostics.get(key) for key in ["n_obs", "n_velocity_valid", "n_sampled", "n_drawn", "reference_magnitude", "display_magnitude_max", "capped_fraction", "input_magnitude_max"]}
        ax.set_title(display + "\n" + "\n".join(labels), fontsize=8, pad=2)
        ax.set_xlabel(""); ax.set_ylabel(""); ax.set_xticks([]); ax.set_yticks([]); ax.axis("off")
    legend_ax = axes.flat[9]
    categories = list(ref.obs["clusters"].cat.categories) if hasattr(ref.obs["clusters"], "cat") else selected
    palette = list(ref.uns.get("clusters_colors", ["#808080"] * len(categories)))
    color_map = dict(zip(categories, palette))
    handles = [Line2D([0], [0], marker="o", color="w", label=cl, markerfacecolor=color_map.get(cl, "#808080"), markersize=7) for cl in selected]
    legend_ax.legend(handles=handles, title="clusters", loc="center", frameon=False, fontsize=9, title_fontsize=10)
    legend_ax.axis("off")
    fig.savefig(out / f"{case_name}_local_phase.png", dpi=300); fig.savefig(out / f"{case_name}_local_phase.pdf"); fig.savefig(out / f"{case_name}_local_phase.svg"); plt.close(fig)
    manifest = {"case": case_name, "gene": case["gene"], "masked_clusters": selected, "edge": case["edge"], "methods": [m[0] for m in METHODS], "edge_scores": score_manifest, "edge_score_source": "evaluation/gene/results/gene_cbdir_edge_details_reference30_all.csv", "layout": "2 rows x 5 columns; 9 methods plus shared cluster legend", "source_paths": [f"artifact://velocity/{p.relative_to(ARTIFACT_ROOT)}" for p in source_paths], "phase_api": "evaluation/gene/plotting/phase_portrait.py:266", "phase_contract": {"layer_s":"Ms", "layer_u":"Mu", "layer_v_s":"velocity", "layer_v_u":"velocity_u", "color_key":"clusters", "grid":[25,25],"density_percentile":15,"arrow_length_mode":"log1p","arrow_reference_quantile":.90,"reference_length":.05,"max_length":.10, "previous_reference_length":.075, "previous_max_length":.15, "s":12, "alpha":0.6, "axes_hidden":True}}
    manifest["phase_diagnostics"] = phase_diagnostics
    manifest["preprocessing_diagnostics"] = preprocessing_diagnostics
    manifest["edge_score_sources"] = ["evaluation/gene/results/gene_cbdir_edge_details_reference30_all.csv"]
    if case["gene"] == "Diaph3" and DIAPH3_SCORES is not None:
        manifest["edge_score_sources"].append("supplement/gene-kinetic/results/diaph3_nipc_to_nbl1.csv")
    (out / f"{case_name}_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(); parser.add_argument("--case", choices=["all", *CASES], default="all"); parser.add_argument("--artifact-root", type=Path, default=None); args = parser.parse_args()
    global ARTIFACT_ROOT
    ARTIFACT_ROOT = args.artifact_root or os.environ.get("VELOCITY_ARTIFACT_ROOT")
    if ARTIFACT_ROOT is None:
        raise SystemExit("Provide --artifact-root or VELOCITY_ARTIFACT_ROOT")
    ARTIFACT_ROOT = Path(ARTIFACT_ROOT).expanduser().resolve()
    sys.path.insert(0, str(REPO))
    from evaluation.gene.run_reference30_gene_cbdir import REAL_CONFIGS
    global REAL_PATHS, DIAPH3_SCORES
    REAL_PATHS = {name: config["paths"] for name, config in REAL_CONFIGS.items()}
    diaph3_path = OUT / "results" / "diaph3_nipc_to_nbl1.csv"
    if diaph3_path.is_file():
        DIAPH3_SCORES = pd.read_csv(diaph3_path)
    from evaluation.gene.plotting import phase_portrait as module
    global EDGE_SCORES
    EDGE_SCORES, _ = load_edge_scores()
    chosen = CASES if args.case == "all" else {args.case: CASES[args.case]}
    for name, case in chosen.items():
        make_panel(name, case, module)


if __name__ == "__main__":
    main()
