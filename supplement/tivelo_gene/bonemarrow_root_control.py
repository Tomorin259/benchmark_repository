#!/usr/bin/env python3
"""Paired post-orientation TIvelo root intervention on archived bone marrow.

Paired sensitivity analysis of TIvelo root selection and trajectory ordering.
Cells, selected genes, moments, undirected neighbors, PAGA and diffusion map
are frozen. Each seed is used for both the archived Mono_1 root and HSC_1.
The canonical two-environment reference-30 cell evaluator scores all outputs.
"""

import argparse
import gc
import hashlib
import importlib.metadata
import json
import random
import subprocess
import sys
import time
from pathlib import Path
from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd
import scvelo as scv
import torch
from scipy.stats import spearmanr
from tivelo.velocity.DTI import get_child_dict, get_d_nn
from tivelo.velocity.model import get_velocity


ROOTS = {"original_root_refit": "Mono_1", "hsc1_root_refit": "HSC_1"}
METRICS = ("cbdir", "batc", "highdim_cbdir", "transcosine", "consistency", "iccoh")
FIT_SETTINGS = dict(
    n_epochs=100, loss_fun="mse", only_s=False, constrain=True,
    alpha_1=1, alpha_2=0.1, batch_size=1024, filter_genes=False,
)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def dense(matrix):
    return matrix.toarray() if hasattr(matrix, "toarray") else np.asarray(matrix)


def children_as_lists(children):
    return {str(key): [str(item) for item in value] for key, value in children.items()}


def save_json(path, payload):
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n")


def check_frozen_inputs(work, original):
    assert work.obs_names.equals(original.obs_names), "Cell order changed"
    assert work.var_names.equals(original.var_names), "Gene order changed"
    for key in ("Mu", "Ms"):
        np.testing.assert_array_equal(dense(work.layers[key]), dense(original.layers[key]))
    for key in ("distances", "connectivities"):
        np.testing.assert_array_equal(dense(work.obsp[key]), dense(original.obsp[key]))
    for key in ("X_diffmap", "X_umap"):
        np.testing.assert_array_equal(work.obsm[key], original.obsm[key])
    np.testing.assert_array_equal(work.uns["diffmap_evals"], original.uns["diffmap_evals"])
    np.testing.assert_array_equal(
        dense(work.uns["paga"]["connectivities"]),
        dense(original.uns["paga"]["connectivities"]),
    )


def graph_audit(work, variant, root, d_nn):
    labels = work.obs["clusters"].astype(str).to_numpy()
    sources, targets = np.nonzero(d_nn)
    cross = labels[sources] != labels[targets]
    counts = pd.DataFrame({"source": labels[sources[cross]], "target": labels[targets[cross]]})
    counts = counts.groupby(["source", "target"]).size().reset_index(name="directed_pairs")
    counts.insert(0, "variant", variant)
    pseudotime = work.obs["dpt_pseudotime"].to_numpy()
    reference = work.obs["palantir_pseudotime"].to_numpy()
    records = []
    for lineage, clusters in (
        ("myeloid", ["HSC_1", "HSC_2", "Precursors", "Mono_1", "Mono_2", "DCs"]),
        ("erythroid", ["HSC_1", "Ery_1", "Ery_2"]),
    ):
        mask = np.isin(labels, clusters) & np.isfinite(pseudotime) & np.isfinite(reference)
        records.append({
            "variant": variant, "root": root, "lineage": lineage,
            "n_cells": int(mask.sum()),
            "spearman_dpt_palantir": float(spearmanr(pseudotime[mask], reference[mask]).statistic),
        })
    return counts, records


def fit_pairs(args, manifest):
    original = ad.read_h5ad(args.source)
    assert str(original.obs.iloc[int(original.uns["iroot"])]["clusters"]) == "Mono_1"
    for key in ("Mu", "Ms", "velocity", "velocity_u"):
        if not np.isfinite(dense(original.layers[key])).all():
            raise ValueError(f"Nonfinite archived layer: {key}")
    output = args.output_dir / "native_outputs"
    logs = args.output_dir / "logs"
    output.mkdir(exist_ok=True)
    logs.mkdir(exist_ok=True)
    torch.set_num_threads(args.threads)
    torch.set_num_interop_threads(1)
    torch.use_deterministic_algorithms(True)
    graph_rows, ordering_rows, templates = [], [], {}
    for variant, root in ROOTS.items():
        work = original.copy()
        # TIvelo recomputes these inside get_child_dict. Freezing them isolates
        # the root intervention and retains the exact archived graph state.
        with patch("scanpy.tl.diffmap", return_value=None), patch("scanpy.tl.paga", return_value=None):
            children, _, _ = get_child_dict(
                work, "clusters", {"start_node": root},
                threshold=0.1, threshold_trans=2.0, adjust=False,
            )
        d_nn, knn, _, same_c = get_d_nn(
            work, "clusters", children, root, root_select="connectivities",
        )
        if variant == "original_root_refit":
            assert children_as_lists(children) == children_as_lists(original.uns["child_dict"])
            np.testing.assert_array_equal(d_nn, dense(original.obsp["d_nn"]))
            np.testing.assert_allclose(
                work.obs["dpt_pseudotime"], original.obs["dpt_pseudotime"],
                atol=1e-8, rtol=1e-8, equal_nan=True,
            )
        check_frozen_inputs(work, original)
        graph, ordering = graph_audit(work, variant, root, d_nn)
        graph_rows.append(graph)
        ordering_rows.extend(ordering)
        manifest["graphs"][variant] = {
            "root": root, "iroot": int(work.uns["iroot"]),
            "root_cell": str(work.obs_names[int(work.uns["iroot"])]),
            "child_dict": children_as_lists(children),
            "directed_neighbor_pairs": int(np.count_nonzero(d_nn)),
            "original_graph_replay_exact": variant == "original_root_refit",
        }
        # The upstream path-selection/orientation step is not rerun here.
        # Avoid carrying an old complete path as if it were newly inferred.
        work.uns["path_dict"] = {"start_node": root}
        for key in list(work.uns):
            if key == "score_main" or key.startswith("score_branch_"):
                del work.uns[key]
        templates[variant] = (work, d_nn, knn, same_c)
    pd.concat(graph_rows, ignore_index=True).to_csv(args.output_dir / "directed_neighbor_pairs.csv", index=False)
    pd.DataFrame(ordering_rows).to_csv(args.output_dir / "dpt_reference_correlations.csv", index=False)
    manifest["n_cells"], manifest["n_genes"] = original.shape
    manifest["archived_path"] = {
        key: (str(value) if isinstance(value, str) else np.asarray(value).tolist())
        for key, value in original.uns["path_dict"].items()
    }
    save_json(args.output_dir / "manifest.json", manifest)

    for seed in args.seeds:
        for variant, root in ROOTS.items():
            stem = f"{variant}_seed{seed}"
            target = output / f"{stem}.h5ad"
            if target.exists():
                if not args.resume or stem not in manifest["fits"]:
                    raise FileExistsError(f"Refusing to overwrite {target}")
                assert sha256(target) == manifest["fits"][stem]["sha256"]
                print(f"SKIP verified {stem}", flush=True)
                continue
            template, d_nn, knn, same_c = templates[variant]
            work = template.copy()
            random.seed(seed)
            np.random.seed(seed)
            torch.manual_seed(seed)
            print(f"FIT START {stem} root={root}", flush=True)
            started = time.monotonic()
            get_velocity(work, d_nn, knn, same_c=same_c, **FIT_SETTINGS)
            for key in ("velocity", "velocity_u"):
                if not np.isfinite(dense(work.layers[key])).all():
                    raise ValueError(f"Nonfinite refitted {key}: {stem}")
            # Never evaluate or save a stale projection from the archived fit.
            for key in ("velocity_graph", "velocity_graph_neg", "velocity_params"):
                work.uns.pop(key, None)
            for key in list(work.obsm):
                if key.startswith("velocity_"):
                    del work.obsm[key]
            scv.tl.velocity_graph(work, vkey="velocity", n_jobs=args.threads)
            scv.tl.velocity_embedding(work, basis="umap", vkey="velocity")
            check_frozen_inputs(work, original)
            work.uns["root_control"] = {
                "variant": variant, "seed": seed, "forced_root": root,
                "intervention_stage": "post-orientation DTI root; frozen PAGA/diffmap",
                "archived_source_sha256": manifest["input_hashes"][str(args.source)],
                "archived_terminal_scores_retained": True,
                "settings_json": json.dumps(FIT_SETTINGS, sort_keys=True),
            }
            work.write_h5ad(target, compression="gzip")
            manifest["fits"][stem] = {
                "variant": variant, "seed": seed, "root": root,
                "path": str(target), "sha256": sha256(target),
                "elapsed_seconds": time.monotonic() - started,
            }
            save_json(args.output_dir / "manifest.json", manifest)
            print(f"FIT DONE {stem} seconds={manifest['fits'][stem]['elapsed_seconds']:.1f}", flush=True)
            del work
            gc.collect()


def evaluate_pairs(args, manifest):
    cell_dir = args.repo_root / "evaluation" / "cell"
    inputs = pd.read_csv(cell_dir / "inputs.tsv", sep="\t", dtype=str, keep_default_na=False)
    bone = inputs[inputs.evaluation_dataset.eq("bonemarrow")]
    native_row = bone[bone.operational_id.eq("TIvelo")]
    reference_row = bone[bone.operational_id.eq("scvelo_deterministic")]
    assert len(native_row) == len(reference_row) == 1
    cases = [{"stem": "archived_original", "variant": "archived_original", "seed": None, "path": str(args.source)}]
    cases.extend({"stem": stem, **record} for stem, record in manifest["fits"].items())
    for case in cases:
        directory = args.output_dir / "evaluations" / case["stem"]
        directory.mkdir(parents=True, exist_ok=True)
        path = Path(case["path"])
        row = native_row.copy()
        row["result_path"] = "velocity:" + path.relative_to(args.artifact_root / "velocity").as_posix()
        row["sha256"] = sha256(path)
        row["file_size_bytes"] = str(path.stat().st_size)
        pd.concat([row, reference_row], ignore_index=True).to_csv(directory / "inputs.tsv", sep="\t", index=False)
        command = [
            sys.executable, str(cell_dir / "run_cell_metrics.py"),
            "--artifact-root", str(args.artifact_root),
            "--inputs", str(directory / "inputs.tsv"),
            "--config", str(cell_dir / "config.json"),
            "--reference-neighbors-root", str(cell_dir / "reference_neighbors"),
            "--reference-moments-root", str(args.reference_moments_root),
            "--output-root", str(directory),
            "--indicators-python", str(args.indicators_python),
            "--tivelo-python", sys.executable,
            "--prepare-script", str(cell_dir / "prepare_cell_metrics.py"),
            "--official-script", str(cell_dir / "score_tivelo_metrics.py"),
            "--method", "TIvelo", "--dataset", "bonemarrow", "--n-jobs", str(args.threads),
        ]
        if args.resume:
            command.append("--resume")
        print(f"EVALUATE START {case['stem']}", flush=True)
        with (args.output_dir / "logs" / f"{case['stem']}.evaluation.log").open("w") as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
        print(f"EVALUATE DONE {case['stem']}", flush=True)
    manifest["evaluation_cases"] = cases
    save_json(args.output_dir / "manifest.json", manifest)


def summarize(args, manifest):
    rows, edges = [], []
    for case in manifest["evaluation_cases"]:
        payload_path = args.output_dir / "evaluations" / case["stem"] / ".work" / "pairs" / "bonemarrow__TIvelo.json"
        payload = json.loads(payload_path.read_text())
        assert payload["status"] == "ok"
        assert payload["metadata"]["result_sha256"] == sha256(Path(case["path"]))
        assert all(item["raw_values_match"] for item in payload["tivelo_parity"])
        rows.append({"variant": case["variant"], "seed": case["seed"], **payload["metrics"]})
        for item in payload["edge_details"] + payload["extended_edge_details"]:
            edges.append({"variant": case["variant"], "seed": case["seed"], **item})
    scores = pd.DataFrame(rows)
    details = pd.DataFrame(edges)
    scores.to_csv(args.output_dir / "cell_scores.csv", index=False)
    details.to_csv(args.output_dir / "edge_scores.csv", index=False)
    canonical = pd.read_csv(args.repo_root / "evaluation/cell/results/cell_metrics.csv")
    supplementary = pd.read_csv(args.repo_root / "evaluation/cell/results/cell_metrics_supplementary.csv")
    canonical = canonical[canonical.dataset.eq("bonemarrow") & canonical.method.eq("TIvelo")].iloc[0]
    supplementary = supplementary[supplementary.dataset.eq("bonemarrow") & supplementary.method.eq("TIvelo")].iloc[0]
    archived = scores[scores.variant.eq("archived_original")].iloc[0]
    differences = {}
    for metric in METRICS:
        reference = supplementary[metric] if metric == "highdim_cbdir" else canonical[metric]
        differences[metric] = abs(float(archived[metric]) - float(reference))
        if differences[metric] > 1e-7:
            raise ValueError(f"Archived reference-30 parity failed: {metric}: {differences[metric]}")
    paired = []
    for seed in args.seeds:
        for metric in METRICS:
            original = float(scores[(scores.seed == seed) & scores.variant.eq("original_root_refit")].iloc[0][metric])
            hsc = float(scores[(scores.seed == seed) & scores.variant.eq("hsc1_root_refit")].iloc[0][metric])
            paired.append({"seed": seed, "metric": metric, "original_root_refit": original,
                           "hsc1_root_refit": hsc, "paired_delta": hsc - original})
    paired = pd.DataFrame(paired)
    paired.to_csv(args.output_dir / "paired_metric_deltas.csv", index=False)
    summary = paired.groupby("metric", sort=False).agg(
        n_pairs=("seed", "size"), original_mean=("original_root_refit", "mean"),
        original_std=("original_root_refit", "std"), hsc1_mean=("hsc1_root_refit", "mean"),
        hsc1_std=("hsc1_root_refit", "std"), delta_mean=("paired_delta", "mean"),
        delta_min=("paired_delta", "min"), delta_max=("paired_delta", "max"),
    ).reset_index()
    summary.to_csv(args.output_dir / "paired_summary.csv", index=False)
    edge_refits = details[details.variant.isin(ROOTS)]
    edge_summary = edge_refits.groupby(["metric", "source", "target", "variant"], sort=False).agg(
        n_seeds=("seed", "size"), mean_score=("score", "mean"), std_score=("score", "std"),
    ).reset_index()
    edge_summary.to_csv(args.output_dir / "edge_summary.csv", index=False)
    manifest["archived_reference30_parity_absolute_differences"] = differences
    manifest["input_hashes_unchanged"] = all(sha256(Path(path)) == expected for path, expected in manifest["input_hashes"].items())
    assert manifest["input_hashes_unchanged"], "A read-only source changed during the experiment"
    manifest["status"] = "complete"
    save_json(args.output_dir / "manifest.json", manifest)
    print(summary.to_string(index=False), flush=True)
    print(edge_summary[edge_summary.metric.eq("cbdir")].to_string(index=False), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-root", type=Path, required=True, help="Parent of the velocity artifact root")
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--reference-moments-root", type=Path, required=True)
    parser.add_argument("--indicators-python", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--seeds", type=int, nargs="+", default=[20261001, 20261002, 20261003])
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--stage", choices=["all", "fit", "evaluate", "summarize"], default="all")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    if len(args.seeds) != len(set(args.seeds)) or len(args.seeds) < 2:
        raise ValueError("Provide at least two distinct seeds for the paired control")
    for name in ("artifact_root", "repo_root", "source", "reference_moments_root", "indicators_python", "output_dir"):
        setattr(args, name, getattr(args, name).resolve())
    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = args.output_dir / "manifest.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        assert manifest["seeds"] == args.seeds
        assert manifest["script_sha256"] == sha256(Path(__file__).resolve())
        for path, expected in manifest["input_hashes"].items():
            assert sha256(Path(path)) == expected, f"Changed input: {path}"
    else:
        cell = args.repo_root / "evaluation/cell"
        frozen = [args.source, cell / "inputs.tsv", cell / "config.json",
                  cell / "reference_neighbors/bonemarrow.npz", args.reference_moments_root / "bonemarrow.npz",
                  cell / "results/cell_metrics.csv", cell / "results/cell_metrics_supplementary.csv",
                  cell / "run_cell_metrics.py", cell / "prepare_cell_metrics.py", cell / "score_tivelo_metrics.py"]
        manifest = {
            "status": "running", "type": "post-orientation root sensitivity analysis",
            "seeds": args.seeds, "fit_settings": FIT_SETTINGS, "threads": args.threads,
            "root_settings": {"t1": 0.1, "t2": 2.0, "root_select": "connectivities (first cell in root cluster)"},
            "held_fixed": ["cells and selected genes", "Mu/Ms", "undirected KNN", "PAGA", "diffusion map", "optimizer and paired seeds", "reference-30 evaluation pipeline"],
            "changed": ["DTI root", "root-dependent child_dict and DPT", "directed neighbors and refitted velocities"],
            "notes": ["Palantir pseudotime is an independent inferred reference, not measured time", "HSC_1 is the biologically specified root in the paired intervention", "X_umap stores the archived t-SNE coordinates in this dataset", "Native velocity projection regenerated before frozen-reference evaluation"],
            "package_versions": {key: importlib.metadata.version(key) for key in ("tivelo", "scvelo", "scanpy", "anndata", "torch", "numpy")},
            "script_sha256": sha256(Path(__file__).resolve()),
            "input_hashes": {str(path): sha256(path) for path in frozen}, "graphs": {}, "fits": {},
        }
        save_json(manifest_path, manifest)
    if args.stage in ("all", "fit"):
        fit_pairs(args, manifest)
    if args.stage in ("all", "evaluate"):
        evaluate_pairs(args, manifest)
    if args.stage in ("all", "summarize"):
        summarize(args, manifest)


if __name__ == "__main__":
    main()
