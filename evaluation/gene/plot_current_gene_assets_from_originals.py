#!/usr/bin/env python3
"""Reuse the established gene-level plotting scripts with current score files."""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import os
import re
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
REPOSITORY_ROOT = HERE.parents[1]
ARTIFACT_ROOT = None
ORIGINAL_REAL = None
ORIGINAL_SIM = None
ORIGINAL_SIM_SUMMARY = None
DEFAULT_OUTPUT = HERE / "figures" / "panels"

sys.path.insert(0, str(REPOSITORY_ROOT))
sys.path.insert(0, str(HERE))
import run_reference30_gene_cbdir as registry  # noqa: E402
from evaluation.gene.plotting import (  # noqa: E402
    PHASE_PLOT_VERSION,
    plot_phase_portrait,
)

DISPLAY_NAMES = {
    "cellDancer": "cellDancer",
    "DeepVelo": "DeepVelo",
    "DeepVelo_scv": "DeepVelo_scv",
    "Graphvelo_scv": "GraphVelo_scv",
    "MultiVelo": "MultiVelo (RNA_only)",
    "scVelo_dyn": "scVelo_dynamic",
    "TIvelo": "TIvelo",
    "UniTVelo": "UniTVelo",
    "veloVI": "veloVI",
}

# Compact labels are used only on dense summary boxplot x-axes. Phase-portrait
# titles and manuscript-facing tables keep the exact naming-contract labels.
SUMMARY_DISPLAY_NAMES = {
    **DISPLAY_NAMES,
    "MultiVelo": "MultiVelo",
    "scVelo_dyn": "scVelo_dyn",
}

# Keep gene-level panels in the same relative order as the final 20-method
# paper order. The remaining methods do not emit both gene-level velocity
# components required by the phase-portrait panels.
GENE_METHOD_ORDER = (
    "scVelo_dyn",
    "MultiVelo",
    "UniTVelo",
    "TIvelo",
    "veloVI",
    "cellDancer",
    "DeepVelo",
    "DeepVelo_scv",
    "Graphvelo_scv",
)
ORIGINAL_REAL_METHOD_ORDER = (
    "cellDancer",
    "DeepVelo",
    "DeepVelo_scv",
    "Graphvelo_scv",
    "MultiVelo",
    "scVelo_dyn",
    "TIvelo",
    "UniTVelo",
    "veloVI",
)
REAL_TXT_KEYS = {
    "cellDancer": "celldancer",
    "DeepVelo": "deepvelo",
    "DeepVelo_scv": "deepvelo2022",
    "Graphvelo_scv": "graphvelo_scv",
    "MultiVelo": "multivelo",
    "scVelo_dyn": "scvelo_dyn",
    "TIvelo": "tivelo",
    "UniTVelo": "unitvelo_with_vu",
    "veloVI": "veloVI",
}
SIM_TXT_KEYS = {
    "cellDancer": "celldancer",
    "DeepVelo": "deepvelo",
    "DeepVelo_scv": "deepvelo2022",
    "Graphvelo_scv": "graphvelo",
    "MultiVelo": "multivelo",
    "scVelo_dyn": "scvelo",
    "TIvelo": "tivelo",
    "UniTVelo": "unitvelo",
    "veloVI": "velovi",
}
REAL_CASE_NAMES = {
    "gastrulation": "Gastrulation",
    "pancreas": "Pancreas",
}
CASE_ORDER = [
    "retina", "DentateGyrus_lamanno", "gastrulation", "capo_pre",
    "bonemarrow", "mousebrain", "pancreas", "branching", "murk", "all",
]
REAL_CASES = [
    "retina", "DentateGyrus_lamanno", "gastrulation", "capo_pre",
    "bonemarrow", "mousebrain", "pancreas",
]
SIM_CASES = ["branching", "murk", "all"]
DATASET_DISPLAY_NAMES = {
    "gastrulation": "Gastrulation",
    "DentateGyrus_lamanno": "Dentategyrus_lamanno",
    "dentategyrus": "Dentategyrus",
    "bonemarrow": "Human bone marrow",
    "capo_pre": "capofibro (VeloCycle)",
    "retina": "Mouse retina",
    "mousebrain": "mousebrain",
    "forebrain": "ForebrainGlut",
    "hematopoiesis": "Human hematopoiesis",
    "pancreas": "Pancreas",
    "scnt_neuron": "scNT_neuron",
    "sceu_organoid": "organoid",
    "sceu_rpe1": "rpe1",
    "rpe1_gen": "rpe1 (VeloCycle)",
    "branching": "Normal + branching",
    "murk": "Normal + MURK",
    "all": "Hybrid (normal + MURK + branching)",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def configure_artifact_root(value: Path) -> None:
    global ARTIFACT_ROOT, ORIGINAL_REAL, ORIGINAL_SIM, ORIGINAL_SIM_SUMMARY
    ARTIFACT_ROOT = Path(value).expanduser().resolve()
    ORIGINAL_REAL = ARTIFACT_ROOT / "plot" / "touse_murk" / "murk_analyze.py"
    ORIGINAL_SIM = ARTIFACT_ROOT / "plot" / "to_use_sim" / "details" / "plot_new_method_comparisons.py"
    ORIGINAL_SIM_SUMMARY = ARTIFACT_ROOT / "plot" / "to_use_sim" / "details" / "plot_new_score_summary.py"


def resolve_artifact_root(value):
    raw_value = value or os.environ.get("VELOCITY_ARTIFACT_ROOT")
    if not raw_value:
        raise SystemExit(
            "Provide --artifact-root or set VELOCITY_ARTIFACT_ROOT to the "
            "external velocity artifact directory."
        )
    return Path(raw_value).expanduser().resolve()


def public_locator(path: Path) -> str:
    path = Path(path)
    if not path.is_absolute():
        return path.as_posix()
    resolved = path.resolve()
    try:
        return resolved.relative_to(REPOSITORY_ROOT.resolve()).as_posix()
    except ValueError:
        pass
    if ARTIFACT_ROOT is not None:
        try:
            relative = resolved.relative_to(ARTIFACT_ROOT.resolve())
        except ValueError:
            relative = None
        if relative is not None:
            return f"artifact://velocity/{relative.as_posix()}"
    return f"artifact://external/{resolved.name}"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def install_repository_phase_plot(module):
    """Keep legacy data loading while routing phase plots through this repo."""
    external_utils = getattr(module, "utils", None)
    if external_utils is None:
        raise AttributeError(f"{module.__name__} does not expose its utils module")
    external_utils.plot_phase_portrait = plot_phase_portrait


def parse_sim_scores(path: Path, case: str) -> pd.DataFrame:
    rows = []
    pattern = re.compile(r"^([^:]+):([^:]+):\s*([^|]+)")
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or "ERROR" in line:
                continue
            match = pattern.match(line.strip())
            cossim_match = re.search(r"CosSim:\s*([^|]+)", line)
            if match is None or cossim_match is None:
                continue
            raw_method, gene, cbdir_text = match.groups()
            method = next(
                (name for name, key in SIM_TXT_KEYS.items() if key == raw_method),
                None,
            )
            if method is not None:
                rows.append({
                    "dataset": case,
                    "method": method,
                    "gene": gene.strip(),
                    "new_cbdir": float(cbdir_text.strip()),
                    "cossim": float(cossim_match.group(1).strip()),
                })
    return pd.DataFrame(rows)


def load_scores(registry, artifact_root):
    root = registry.HERE.parent / "results"
    score_path = root / "gene_cbdir_reference30_all.csv"
    score_paths = {case: score_path for case in CASE_ORDER}
    frames = []
    for case in CASE_ORDER:
        frame = pd.read_csv(score_paths[case])
        frames.append(frame[frame["dataset"].eq(case)][["dataset", "method", "gene", "new_cbdir"]])
    scores = pd.concat(frames, ignore_index=True)
    compact_path = HERE / "results" / "simulation" / "gene_cbdir_cossim_reference.csv"
    if compact_path.is_file():
        compact = pd.read_csv(compact_path)
        required = {"Dataset", "Gene", "Method", "CBDir", "CosSim"}
        missing = required.difference(compact.columns)
        if missing:
            raise ValueError(f"Compact CosSim table is missing columns: {sorted(missing)}")
        sim = compact.rename(columns={
            "Dataset": "dataset", "Gene": "gene", "Method": "method",
            "CBDir": "new_cbdir", "CosSim": "cossim",
        })[["dataset", "method", "gene", "new_cbdir", "cossim"]]
        sim_paths = {case: compact_path for case in SIM_CASES}
    else:
        sim_paths = {
            "branching": artifact_root / "indicators/gene/simulation/branching_analyze/results_new/score_rerun_results.txt",
            "murk": artifact_root / "indicators/gene/simulation/murk_analyze/results_new/score_rerun_results.txt",
            "all": artifact_root / "indicators/gene/simulation/sim_all_analyze/results_new/score_rerun_results.txt",
        }
        sim = pd.concat([parse_sim_scores(path, case) for case, path in sim_paths.items()], ignore_index=True)
    scores = scores.merge(sim[["dataset", "method", "gene", "cossim"]], on=["dataset", "method", "gene"], how="left")
    return scores, {**score_paths, **sim_paths}


def configure_real_module(module):
    module.METHOD_NAMES = [DISPLAY_NAMES[name] for name in GENE_METHOD_ORDER]
    module.TXT_METHOD_KEYS = [REAL_TXT_KEYS[name] for name in GENE_METHOD_ORDER]
    old_colors = module.BENCHMARK_COLORS
    module.BENCHMARK_COLORS = {
        DISPLAY_NAMES.get(old_name, old_name): color
        for old_name, color in old_colors.items()
    }
    return list(GENE_METHOD_ORDER)


def write_real_score_file(path: Path, case: str, config, scores: pd.DataFrame, txt_keys):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for method, txt_key in zip(GENE_METHOD_ORDER, txt_keys):
            for gene in config.genes:
                row = scores[(scores.dataset == case) & (scores.method == method) & (scores.gene == gene)]
                if row.empty or not np.isfinite(row.iloc[0]["new_cbdir"]):
                    continue
                handle.write(
                    f"{txt_key}:{gene}: {row.iloc[0]['new_cbdir']:.5f} | Raw_Scores: {{}}\n"
                )


def write_sim_score_file(path: Path, case: str, scores: pd.DataFrame):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for method in GENE_METHOD_ORDER:
            raw_key = SIM_TXT_KEYS[method]
            for gene in registry.SIM_CONFIGS[case]["genes"]:
                row = scores[(scores.dataset == case) & (scores.method == method) & (scores.gene == gene)]
                if row.empty or not np.isfinite(row.iloc[0]["new_cbdir"]):
                    continue
                cossim = row.iloc[0]["cossim"]
                handle.write(
                    f"{raw_key}:{gene}: {row.iloc[0]['new_cbdir']:.5f} | "
                    f"CosSim: {cossim:.5f} | Raw_Scores: {{}}\n"
                )


def write_scores_used(path: Path, case: str, scores: pd.DataFrame):
    frame = scores[scores.dataset == case].copy()
    frame["display_method"] = frame["method"].map(DISPLAY_NAMES)
    frame["display_dataset"] = DATASET_DISPLAY_NAMES.get(case, case)
    frame.to_csv(path, index=False)


def set_real_summary_style(real_module):
    saved_names = real_module.METHOD_NAMES
    saved_colors = real_module.BENCHMARK_COLORS
    real_module.METHOD_NAMES = [SUMMARY_DISPLAY_NAMES[name] for name in GENE_METHOD_ORDER]
    real_module.BENCHMARK_COLORS = {
        SUMMARY_DISPLAY_NAMES[name]: saved_colors.get(DISPLAY_NAMES[name], "#808080")
        for name in GENE_METHOD_ORDER
    }
    return saved_names, saved_colors


def render_real_case(real_module, case: str, scores: pd.DataFrame, output_root: Path):
    config_name = REAL_CASE_NAMES.get(case, case)
    source_config = real_module.DATASET_CONFIGS[config_name]
    path_by_method = dict(zip(ORIGINAL_REAL_METHOD_ORDER, source_config.h5ad_paths))
    current_paths = [
        path_by_method[method].replace(
            "velocity/results/TIvelo/TIvelo_for_cells/", "velocity/results/TIvelo/"
        )
        for method in GENE_METHOD_ORDER
    ]
    output_dir = output_root / "real" / case
    config = replace(source_config, h5ad_paths=current_paths, output_dir=str(output_dir))
    write_real_score_file(config.score_txt_path, case, config, scores, real_module.TXT_METHOD_KEYS)
    adatas = real_module.load_and_preprocess_dataset(config)
    real_module.plot_multi_method_genes(adatas, config)
    saved_names, saved_colors = set_real_summary_style(real_module)
    real_module.plot_score_summary_box_heatmap(config)
    real_module.METHOD_NAMES = saved_names
    real_module.BENCHMARK_COLORS = saved_colors
    summary = output_dir / "branching_combined.png"
    if summary.is_file():
        summary.rename(output_dir / "summary_cbdir.png")
    write_scores_used(output_dir / "scores_used.csv", case, scores)
    del adatas


def render_sim_cases(sim_module, summary_module, cases, scores, output_root):
    sim_root = output_root / "simulation"
    score_root = output_root / "score_sources"
    sim_root.mkdir(parents=True, exist_ok=True)
    score_root.mkdir(parents=True, exist_ok=True)
    sim_module.METHODS = [
        (SIM_TXT_KEYS[name], DISPLAY_NAMES[name]) for name in GENE_METHOD_ORDER
    ]
    sim_summary_labels = {SIM_TXT_KEYS[name]: SUMMARY_DISPLAY_NAMES[name] for name in GENE_METHOD_ORDER}
    sim_summary_module = summary_module
    sim_summary_module.METHOD_LABELS = sim_summary_labels
    sim_summary_module.PALETTE = {SUMMARY_DISPLAY_NAMES[name]: color for name, color in {
        "cellDancer": "#5BA658", "DeepVelo": "#C9383B", "DeepVelo_scv": "#5F6B6D",
        "Graphvelo_scv": "#B5BD61", "MultiVelo": "#8C564B", "scVelo_dyn": "#9467BD",
        "TIvelo": "#D16FA3", "UniTVelo": "#E377C2", "veloVI": "#BCBD22",
    }.items()}
    for case in cases:
        score_path = score_root / f"{case}.txt"
        write_sim_score_file(score_path, case, scores)
        sim_module.SCORE_PATHS[case] = score_path
        sim_module.PLOT_BASE = sim_root
        sim_module.plot_dataset(case)

        summary_module.SCORE_FILES[case] = score_path
        summary_df = pd.DataFrame(summary_module.parse_score_file(case, score_path))
        summary_df["display_method"] = summary_df["Method"]
        summary_df["display_dataset"] = DATASET_DISPLAY_NAMES.get(case, case)
        case_dir = sim_root / case
        summary_df.to_csv(case_dir / "scores_used.csv", index=False)
        summary_module.plot_boxplots(summary_df, case_dir / "summary_cbdir_cossim.png")


def render_summaries_only(real_module, summary_module, real_cases, sim_cases, output_root: Path):
    """Refresh summary boxplots while preserving existing phase panels."""
    for case in real_cases:
        config_name = REAL_CASE_NAMES.get(case, case)
        source_config = real_module.DATASET_CONFIGS[config_name]
        output_dir = output_root / "real" / case
        config = replace(source_config, output_dir=str(output_dir))
        saved_names, saved_colors = set_real_summary_style(real_module)
        real_module.plot_score_summary_box_heatmap(config)
        real_module.METHOD_NAMES = saved_names
        real_module.BENCHMARK_COLORS = saved_colors
        summary = output_dir / "branching_combined.png"
        if summary.is_file():
            summary.replace(output_dir / "summary_cbdir.png")

    sim_root = output_root / "simulation"
    score_root = output_root / "score_sources"
    sim_summary_module = summary_module
    sim_summary_module.METHOD_LABELS = {
        SIM_TXT_KEYS[name]: SUMMARY_DISPLAY_NAMES[name] for name in GENE_METHOD_ORDER
    }
    sim_summary_module.PALETTE = {SUMMARY_DISPLAY_NAMES[name]: color for name, color in {
        "cellDancer": "#5BA658", "DeepVelo": "#C9383B", "DeepVelo_scv": "#5F6B6D",
        "Graphvelo_scv": "#B5BD61", "MultiVelo": "#8C564B", "scVelo_dyn": "#9467BD",
        "TIvelo": "#D16FA3", "UniTVelo": "#E377C2", "veloVI": "#BCBD22",
    }.items()}
    for case in sim_cases:
        score_path = score_root / f"{case}.txt"
        sim_summary_module.SCORE_FILES[case] = score_path
        summary_df = pd.DataFrame(sim_summary_module.parse_score_file(case, score_path))
        sim_summary_module.plot_boxplots(summary_df, sim_root / case / "summary_cbdir_cossim.png")


def write_manifest(output_root: Path, score_sources: dict, cases):
    rows = []
    mapping = {
        "Figure4_B": [("murk", gene) for gene in registry.SIM_CONFIGS["murk"]["genes"]],
        "Figure4_C": [("branching", gene) for gene in registry.SIM_CONFIGS["branching"]["genes"]],
        "Figure4_D": [(case, "summary") for case in SIM_CASES],
        "Figure5_A": [("gastrulation", "gene_phase_raw")],
        "Figure5_B": [("gastrulation", gene) for gene in registry.REAL_CONFIGS["gastrulation"]["genes"]],
        "Figure5_C": [("gastrulation", "summary")],
        "Figure5_E": [(case, "gene_phase_raw") for case in ("DentateGyrus_lamanno", "retina")],
        "Figure5_G": [("retina", gene) for gene in registry.REAL_CONFIGS["retina"]["genes"]],
        "Figure5_F": [("DentateGyrus_lamanno", gene) for gene in registry.REAL_CONFIGS["DentateGyrus_lamanno"]["genes"]],
        "Figure5_H": [(case, "summary") for case in ("retina", "DentateGyrus_lamanno")],
        "SuppFig1": [("DentateGyrus_lamanno", gene) for gene in registry.REAL_CONFIGS["DentateGyrus_lamanno"]["genes"] if gene != "Ncald"],
        "SuppFig2": [("retina", gene) for gene in registry.REAL_CONFIGS["retina"]["genes"] if gene != "Snap25"],
        "SuppFig3": [("gastrulation", gene) for gene in registry.REAL_CONFIGS["gastrulation"]["genes"][1:11]],
        "SuppFig4": [("gastrulation", gene) for gene in registry.REAL_CONFIGS["gastrulation"]["genes"][11:]],
        "SuppFig5": [("capo_pre", gene) for gene in registry.REAL_CONFIGS["capo_pre"]["genes"]],
        "SuppFig6": [("bonemarrow", gene) for gene in registry.REAL_CONFIGS["bonemarrow"]["genes"]],
        "SuppFig7": [("mousebrain", gene) for gene in registry.REAL_CONFIGS["mousebrain"]["genes"]],
        "SuppFig8": [("pancreas", gene) for gene in registry.REAL_CONFIGS["pancreas"]["genes"]],
    }
    for panel, assets in mapping.items():
        for case, asset in assets:
            root = output_root / ("simulation" if case in SIM_CASES else "real") / case
            if asset == "summary":
                candidates = list(root.glob("summary*.png"))
            elif asset == "gene_phase_raw":
                candidates = [root / "gene_phase_raw.png"]
            else:
                candidates = [root / f"{asset}_9methods_compare.png"]
            for path in candidates:
                if not path.is_file():
                    continue
                rows.append({
                    "panel": panel, "dataset": case, "asset": asset,
                    "path": str(path.relative_to(output_root)),
                    "score_source": (
                        public_locator(score_sources[case])
                        if asset != "gene_phase_raw"
                        else ""
                    ),
                    "sha256": sha256(path),
                    "source_kind": (
                        "phase_portrait_panel"
                        if asset == "gene_phase_raw"
                        else "score_annotated_panel"
                    ),
                })
    with (output_root / "assembly_manifest.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--artifact-root",
        type=Path,
        default=None,
        help="External velocity artifact directory (or VELOCITY_ARTIFACT_ROOT).",
    )
    parser.add_argument("--cases", nargs="+", choices=CASE_ORDER, default=CASE_ORDER)
    parser.add_argument(
        "--summaries-only",
        action="store_true",
        help="Refresh only summary boxplots; keep existing phase-portrait panels.",
    )
    args = parser.parse_args()
    configure_artifact_root(resolve_artifact_root(args.artifact_root))
    registry.ARTIFACT_ROOT = ARTIFACT_ROOT
    registry.REFERENCE_NEIGHBOR_ROOT = ARTIFACT_ROOT / "indicators" / "cell" / "reference_neighbors"
    output_root = (
        args.output_root
        if args.output_root.is_absolute()
        else REPOSITORY_ROOT / args.output_root
    ).resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    os.chdir(ARTIFACT_ROOT.parent)

    scores, score_sources = load_scores(registry, ARTIFACT_ROOT)
    sys.path.insert(0, str(ORIGINAL_REAL.parent))
    real_module = load_module(ORIGINAL_REAL, "murk_analyze_original_for_current_assets")
    install_repository_phase_plot(real_module)
    configure_real_module(real_module)
    sim_plot_dir = ORIGINAL_SIM.parent
    sys.path.insert(0, str(sim_plot_dir))
    sys.modules.pop("utils", None)
    sim_module = load_module(ORIGINAL_SIM, "simulation_plot_original_for_current_assets")
    install_repository_phase_plot(sim_module)
    sys.modules.pop("utils", None)
    summary_module = load_module(ORIGINAL_SIM_SUMMARY, "simulation_summary_original_for_current_assets")

    selected_real = [case for case in args.cases if case in REAL_CASES]
    selected_sim = [case for case in args.cases if case in SIM_CASES]
    if args.summaries_only:
        render_summaries_only(real_module, summary_module, selected_real, selected_sim, output_root)
    else:
        for case in selected_real:
            render_real_case(real_module, case, scores, output_root)
        if selected_sim:
            render_sim_cases(sim_module, summary_module, selected_sim, scores, output_root)

    provenance = {
        "artifact_root": "artifact://velocity",
        "renderer": "evaluation/gene/plot_current_gene_assets_from_originals.py",
        "renderer_sha256": sha256(HERE / "plot_current_gene_assets_from_originals.py"),
        "original_real_plot_script": "artifact://velocity/plot/touse_murk/murk_analyze.py",
        "original_real_plot_script_sha256": sha256(ORIGINAL_REAL),
        "original_sim_plot_script": "artifact://velocity/plot/to_use_sim/details/plot_new_method_comparisons.py",
        "original_sim_plot_script_sha256": sha256(ORIGINAL_SIM),
        "original_sim_summary_script": "artifact://velocity/plot/to_use_sim/details/plot_new_score_summary.py",
        "original_sim_summary_script_sha256": sha256(ORIGINAL_SIM_SUMMARY),
        "cases": args.cases,
        "display_method_names": DISPLAY_NAMES,
        "summary_display_method_names": SUMMARY_DISPLAY_NAMES,
        "display_dataset_names": DATASET_DISPLAY_NAMES,
        "phase_plot_renderer": "evaluation/gene/plotting/phase_portrait.py",
        "phase_plot_renderer_sha256": sha256(REPOSITORY_ROOT / "evaluation/gene/plotting/phase_portrait.py"),
        "phase_plot_contract": {
            "version": PHASE_PLOT_VERSION,
            "dpi": 600,
            "arrow_grid": [25, 25],
            "layer_s": "Ms",
            "layer_u": "Mu",
            "arrow_length_mode": "log1p",
            "arrow_reference_quantile": 0.90,
            "arrow_reference_length": 0.05,
            "arrow_max_length": 0.10,
            "arrow_density_percentile": 15.0,
            "randomized_sampling": False,
            "velocity_smoothing": False,
        },
        "score_sources": {
            key: public_locator(value) for key, value in score_sources.items()
        },
    }
    (output_root / "provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_manifest(output_root, score_sources, args.cases)
    print(f"Wrote original-script gene figure assets to {output_root}")


if __name__ == "__main__":
    main()
