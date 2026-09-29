#!/usr/bin/env python3
"""Run the official TSvelo workflow on the 14 benchmark datasets."""

from pathlib import Path
import argparse
import os
import shutil
import subprocess
import sys
import tempfile


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = REPOSITORY_ROOT.parent / "data"
OUTPUT_ROOT = REPOSITORY_ROOT / "outputs" / "tsvelo"

DATASETS_TO_RUN = [
    "bonemarrow", "capo_pre", "DentateGyrus_lamanno", "ForebrainGlut",
    "Gastrulation", "Pancreas", "dentategyrus", "hematopoiesis_dynamo",
    "mousebrain", "neuron", "organoid", "retina", "rpe1", "rpe1_gen",
]

# input, label, official dataset name, cell-ID normalization
DATASETS = {
    "bonemarrow": ("bonemarrow/bonemarrow.h5ad", "clusters", "bonemarrow", False),
    "capo_pre": ("capo_pre/capo_pre.h5ad", "phase", "capofibro", False),
    "DentateGyrus_lamanno": (
        "DentateGyrus_lamanno/DentateGyrus_lamanno.h5ad",
        "clusters",
        "DentateGyrus_lamanno_1",
        False,
    ),
    "ForebrainGlut": ("ForebrainGlut/ForebrainGlut.loom", "Clusters", "forebrain", False),
    "Gastrulation": ("Gastrulation/Gastrulation.h5ad", "celltype", "Gastrulation_ori", False),
    "Pancreas": ("Pancreas/Pancreas.h5ad", "clusters", "Pancreas", False),
    "dentategyrus": ("dentategyrus/dentategyrus.h5ad", "clusters", "dentategyrus", False),
    "hematopoiesis_dynamo": (
        "hematopoiesis_dynamo/hematopoiesis_dynamo.h5ad",
        "cell_type",
        "hematopoiesis",
        False,
    ),
    "mousebrain": (
        "mousebrain/10X_multiome_mouse_brain.loom",
        "celltype",
        "mousebrain",
        True,
    ),
    "neuron": ("neuron/neuron.h5ad", "time_str", "neuron", False),
    "organoid": ("organoid/organoid_pre.h5ad", "cell_type", "organoid_ori", False),
    "retina": ("retina/retina_full.h5ad", "Annotation", "retina_ori", False),
    "rpe1": ("rpe1/rpe1.h5ad", "cell_cycle_phase", "rpe1", False),
    "rpe1_gen": ("rpe1_gen/rpe1_gen.h5ad", "phase", "rpe1_gen", False),
}


def get_source(path):
    path = path or os.environ.get("TSVELO_ROOT")
    if not path:
        raise ValueError("set --source-root or TSVELO_ROOT")
    return Path(path).expanduser().resolve()


def check_paths(source, data_root, datasets):
    required = (
        source / "TSvelo_run.py",
        source / "TSvelo" / "TSvelo_pp.py",
        source / "TSvelo" / "TSvelo_branch.py",
        source / "TSvelo" / "TSvelo_concate.py",
        source / "ENCODE" / "processed",
        source / "ChEA" / "ChEA_2016.txt",
    )
    missing = [str(path) for path in required if not path.exists()]
    missing += [
        str(data_root / DATASETS[dataset][0])
        for dataset in datasets
        if not (data_root / DATASETS[dataset][0]).exists()
    ]
    if missing:
        raise FileNotFoundError("\n".join(missing))


def normalize_cell_id(value):
    value = str(value)
    value = value.rsplit(":", 1)[-1]
    return value[:-1] + "-1" if value.endswith("x") else value


def stage_input(input_path, label_key, target_path, normalize_ids=False):
    import scanpy as sc

    adata = sc.read_loom(input_path, sparse=True) if input_path.suffix == ".loom" else sc.read_h5ad(input_path)
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    if label_key not in adata.obs:
        raise KeyError(f"{input_path}: missing obs[{label_key!r}]")
    # The official fallback loader copies cell_type into clusters.
    adata.obs["cell_type"] = adata.obs[label_key].copy()
    if normalize_ids:
        adata.obs_names = [normalize_cell_id(value) for value in adata.obs_names]
        adata.obs_names_make_unique()
    target_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(target_path)


def run_dataset(dataset, source, data_root, output_root, cuda=2, n_jobs=-1, overwrite=False):
    input_path, label_key, official_name, normalize_ids = DATASETS[dataset]
    destination = output_root / official_name
    if destination.exists():
        if not overwrite:
            raise FileExistsError(f"output exists: {destination}")
        shutil.rmtree(destination)

    with tempfile.TemporaryDirectory(prefix=f"tsvelo_{dataset}_") as temp:
        work = Path(temp)
        package = work / "velocity" / "methods" / "TSvelo"
        package.mkdir(parents=True)
        (package / "ENCODE").symlink_to(source / "ENCODE", target_is_directory=True)
        (package / "ChEA").symlink_to(source / "ChEA", target_is_directory=True)
        (package / "results").mkdir()
        stage_input(
            data_root / input_path,
            label_key,
            work / "velocity" / "data" / "organoid" / "organoids_processed_for_tfvelo.h5ad",
            normalize_ids,
        )
        env = os.environ.copy()
        env["PYTHONPATH"] = os.pathsep.join([str(source), env.get("PYTHONPATH", "")]).rstrip(os.pathsep)
        command = [
            sys.executable,
            str(source / "TSvelo_run.py"),
            "--dataset_name", official_name,
            "--n_jobs", str(n_jobs),
            "--n_neighbors", "30",
            "--n_top_genes", "2000",
            "--n_selected_genes", "100",
            "--TF_databases", "ENCODE", "ChEA",
            "--N_steps", "1000",
            "--cuda", str(cuda),
            "--N_EPOCH", "30",
            "--num_epochs", "500",
            "--min_decrease", "0",
            "--n_genes2show", "0",
        ]
        subprocess.run(command, cwd=work, env=env, check=True)
        result = package / "results" / official_name
        if not (result / "TSvelo.h5ad").exists():
            raise RuntimeError(f"{dataset}: TSvelo.h5ad was not created")
        output_root.mkdir(parents=True, exist_ok=True)
        shutil.copytree(result, destination)
    print(f"wrote: {destination}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", choices=[*DATASETS_TO_RUN, "all"], default="all")
    parser.add_argument("--source-root", help="official TSvelo checkout; defaults to TSVELO_ROOT")
    parser.add_argument("--data-root", default=str(DATA_ROOT))
    parser.add_argument("--output-root", default=str(OUTPUT_ROOT))
    parser.add_argument("--cuda", type=int, default=2)
    parser.add_argument("--n-jobs", type=int, default=-1)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()

    source = get_source(args.source_root)
    data_root = Path(args.data_root).expanduser().resolve()
    datasets = DATASETS_TO_RUN if args.dataset == "all" else [args.dataset]
    check_paths(source, data_root, datasets)
    if args.check_only:
        print(f"TSvelo configuration valid: {len(datasets)} datasets", flush=True)
        return
    for dataset in datasets:
        run_dataset(
            dataset,
            source,
            data_root,
            Path(args.output_root).expanduser().resolve(),
            cuda=args.cuda,
            n_jobs=args.n_jobs,
            overwrite=args.overwrite,
        )


if __name__ == "__main__":
    main()
