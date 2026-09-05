#!/usr/bin/env python3
"""Run the official TFvelo demo for the benchmark datasets."""

from pathlib import Path
import argparse
import os
import shutil
import subprocess
import sys
import tempfile


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = REPOSITORY_ROOT.parent / "data"
OUTPUT_ROOT = REPOSITORY_ROOT / "outputs" / "tfvelo"

DATASETS_TO_RUN = [
    "bonemarrow", "capo_pre", "DentateGyrus_lamanno", "ForebrainGlut",
    "Gastrulation", "Pancreas", "dentategyrus", "hematopoiesis_dynamo",
    "mousebrain", "neuron", "organoid", "retina", "rpe1", "rpe1_gen",
]

# input, label, official dataset name, plotting basis
DATASETS = {
    "bonemarrow": ("bonemarrow/bonemarrow.h5ad", "clusters", "bonemarrow", "umap"),
    "capo_pre": ("capo_pre/capo_pre.h5ad", "phase", "CapoFibro", "umap"),
    "DentateGyrus_lamanno": (
        "DentateGyrus_lamanno/DentateGyrus_lamanno.h5ad", "clusters",
        "DentateGyrus_lamanno", "umap",
    ),
    "ForebrainGlut": ("ForebrainGlut/ForebrainGlut.loom", "Clusters", "forebrain", "pca"),
    "Gastrulation": ("Gastrulation/Gastrulation.h5ad", "celltype", "gastrulation_erythroid", "umap"),
    "Pancreas": ("Pancreas/Pancreas.h5ad", "clusters", "pancreas", "umap"),
    "dentategyrus": ("dentategyrus/dentategyrus.h5ad", "clusters", "dentategyrus", "umap"),
    "hematopoiesis_dynamo": (
        "hematopoiesis_dynamo/hematopoiesis_dynamo.h5ad", "cell_type",
        "hematopoiesis", "umap",
    ),
    "mousebrain": ("mousebrain/mousebrain.h5ad", "celltype", "mousebrain", "umap"),
    "neuron": ("neuron/neuron.h5ad", "time_str", "neuron", "umap"),
    "organoid": ("organoid/organoid_pre.h5ad", "cell_type", "organoid", "umap"),
    "retina": ("retina/retina_full.h5ad", "Annotation", "retina", "umap"),
    "rpe1": ("rpe1/rpe1.h5ad", "cell_cycle_phase", "rpe1", "umap"),
    "rpe1_gen": ("rpe1_gen/rpe1_gen.h5ad", "phase", "rpe1_gen", "umap"),
}

STAGE_PATHS = {
    "Gastrulation": "velocity/data/Gastrulation/erythroid_lineage.h5ad",
    "Pancreas": "velocity/data/Pancreas/endocrinogenesis_day15.h5ad",
}


def get_source(path):
    path = path or os.environ.get("TFVELO_ROOT")
    if not path:
        raise ValueError("set --source-root or TFVELO_ROOT")
    return Path(path).expanduser().resolve()


def check_paths(source, data_root):
    required = (
        source / "TFvelo_run_demo.py",
        source / "TFvelo_analysis_demo.py",
        source / "TFvelo" / "__init__.py",
        source / "ENCODE" / "processed",
        source / "ChEA" / "ChEA_2016.txt",
    )
    missing = [str(path) for path in required if not path.exists()]
    missing += [
        str(data_root / DATASETS[dataset][0])
        for dataset in DATASETS_TO_RUN
        if not (data_root / DATASETS[dataset][0]).exists()
    ]
    if missing:
        raise FileNotFoundError("\n".join(missing))


def stage_input(input_path, label_key, target_path):
    import scanpy as sc

    adata = sc.read_loom(input_path, sparse=True) if input_path.suffix == ".loom" else sc.read_h5ad(input_path)
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    if label_key not in adata.obs:
        raise KeyError(f"{input_path}: missing obs[{label_key!r}]")
    adata.obs["clusters"] = adata.obs[label_key].copy()
    target_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(target_path)


def patch_analysis(source, target):
    old = 'adata = ad.read_h5ad("TFvelo_results/TFvelo_hematopoiesis_demo/rc.h5ad")'
    text = (source / "TFvelo_analysis_demo.py").read_text(encoding="utf-8")
    if text.count(old) != 1:
        raise RuntimeError("TFvelo_analysis_demo.py input line changed upstream")
    target.write_text(
        text.replace(old, "adata = ad.read_h5ad(args.data_path + 'rc.h5ad')"),
        encoding="utf-8",
    )


def run_dataset(dataset, source, data_root, output_root, overwrite=False):
    input_path, label_key, official_name, basis = DATASETS[dataset]
    destination = output_root / f"TFvelo_{official_name}_demo"
    if destination.exists():
        if not overwrite:
            raise FileExistsError(f"output exists: {destination}")
        shutil.rmtree(destination)

    with tempfile.TemporaryDirectory(prefix=f"tfvelo_{dataset}_") as temp:
        work = Path(temp)
        (work / "TFvelo-main").symlink_to(source, target_is_directory=True)
        (work / "velocity" / "data").mkdir(parents=True)
        staged = work / "velocity" / "data" / "dentategyrus" / "dentategyrus.h5ad"
        if dataset in STAGE_PATHS:
            staged = work / STAGE_PATHS[dataset]
        stage_input(data_root / input_path, label_key, staged)
        (work / "TFvelo_results").mkdir()
        analysis = work / "TFvelo_analysis_demo_for_benchmark.py"
        patch_analysis(source, analysis)
        env = os.environ.copy()
        env["PYTHONPATH"] = os.pathsep.join([str(source), env.get("PYTHONPATH", "")]).rstrip(os.pathsep)
        args = ["--dataset_name", official_name, "--basis", basis]
        subprocess.run([sys.executable, str(source / "TFvelo_run_demo.py"), *args], cwd=work, env=env, check=True)
        subprocess.run([sys.executable, str(analysis), *args], cwd=work, env=env, check=True)
        result = work / "TFvelo_results" / destination.name
        if not (result / "TFvelo.h5ad").exists():
            raise RuntimeError(f"{dataset}: TFvelo.h5ad was not created")
        output_root.mkdir(parents=True, exist_ok=True)
        shutil.copytree(result, destination)
    print(f"wrote: {destination}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", choices=[*DATASETS_TO_RUN, "all"], default="all")
    parser.add_argument("--source-root", help="official TFvelo checkout; defaults to TFVELO_ROOT")
    parser.add_argument("--data-root", default=str(DATA_ROOT))
    parser.add_argument("--output-root", default=str(OUTPUT_ROOT))
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()

    source = get_source(args.source_root)
    data_root = Path(args.data_root).expanduser().resolve()
    check_paths(source, data_root)
    if args.check_only:
        print(f"TFvelo configuration valid: {len(DATASETS)} datasets", flush=True)
        return
    datasets = DATASETS_TO_RUN if args.dataset == "all" else [args.dataset]
    for dataset in datasets:
        run_dataset(dataset, source, data_root, Path(args.output_root).expanduser().resolve(), args.overwrite)


if __name__ == "__main__":
    main()
