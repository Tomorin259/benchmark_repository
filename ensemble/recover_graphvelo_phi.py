#!/usr/bin/env python3
"""Verify saved or replay native GraphVelo-SCV coefficients from registered inputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
from importlib.metadata import version
import json
from pathlib import Path

import anndata as ad
from anndata._io.specs import read_elem
from graphvelo.graph_velocity import (
    GraphVelo, corr_kernel, cos_corr, density_corrected_transition_matrix,
    tangent_space_projection,
)
import h5py
import numpy as np
from scipy import sparse
import graphvelo.graph_velocity as graph_velocity
from threadpoolctl import threadpool_info


ROOT = Path(__file__).resolve().parents[1]
ATOL = 1e-7
RTOL = 1e-6
ZERO_TOLERANCE = 1e-12


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def compare_arrays(actual, expected):
    actual = actual.toarray() if sparse.issparse(actual) else np.asarray(actual)
    expected = expected.toarray() if sparse.issparse(expected) else np.asarray(expected)
    if actual.shape != expected.shape:
        raise ValueError("replayed and stored output shapes differ")
    actual = actual.astype(np.float64, copy=False)
    expected = expected.astype(np.float64, copy=False)
    finite = np.isfinite(actual) & np.isfinite(expected)
    equal_invalid = np.array_equal(np.isfinite(actual), np.isfinite(expected))
    delta = actual[finite] - expected[finite]
    reference_norm = np.linalg.norm(expected[finite])
    return {
        "shape": list(actual.shape),
        "max_abs_error": float(np.max(np.abs(delta))) if delta.size else 0.0,
        "relative_l2_error": float(np.linalg.norm(delta) / max(reference_norm, 1e-30)),
        "nonfinite_pattern_matches": bool(equal_invalid),
        "entries_outside_tolerance": int(np.count_nonzero(
            np.abs(delta) > ATOL + RTOL * np.abs(expected[finite])
        )),
        "allclose": bool(np.allclose(actual, expected, atol=ATOL, rtol=RTOL, equal_nan=True)),
        "entrywise_equal": bool(np.array_equal(actual, expected, equal_nan=True)),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--inputs", type=Path, default=ROOT / "evaluation/cell/inputs.tsv")
    parser.add_argument("--dataset", action="append")
    parser.add_argument("--n-jobs", type=int, default=4)
    parser.add_argument("--mack-top", type=int, help="explicit historical MacK-ranked gene count")
    args = parser.parse_args()
    if args.n_jobs < 1 or (args.mack_top is not None and args.mack_top < 1):
        raise ValueError("--n-jobs and --mack-top must be positive")
    with args.inputs.open(newline="") as handle:
        registered = [r for r in csv.DictReader(handle, delimiter="\t")
                      if r["operational_id"] == "graphvelo_scv"]
    selected = args.dataset or [r["evaluation_dataset"] for r in registered]
    if len(set(selected)) != len(selected) or set(selected).difference(
        r["evaluation_dataset"] for r in registered
    ):
        raise ValueError("select unique registered GraphVelo-SCV dataset keys")
    args.output_root.mkdir(parents=True, exist_ok=True)
    source_file = Path(graph_velocity.__file__)
    protocol = {
        "graphvelo_version": version("GraphVelo"),
        "graphvelo_source_sha256": sha256(source_file),
        "generator_sha256": sha256(Path(__file__)),
        "input_registry_sha256": sha256(args.inputs),
        "initialization": {"xkey": "Ms", "vkey": "velocity",
                           "gene_subset": "mack_top_" + str(args.mack_top) if args.mack_top else None,
                           "approx": True, "n_pcs": 30},
        "fit": {"a": 1, "b": 10, "r": 1, "loss_func": "linear",
                "softmax_adjusted": False, "n_jobs": args.n_jobs},
        "verification": {"atol": ATOL, "rtol": RTOL},
        "normalization": {
            "nonzero": "offdiagonal Phi divided by its row signed-L1 norm; diagonal zero",
            "zero": "identity row if the offdiagonal signed-L1 norm is at most 1e-12",
        },
        "threadpools": [{key: pool.get(key) for key in
                         ("user_api", "internal_api", "num_threads", "version", "architecture")}
                        for pool in threadpool_info()],
        "versions": {name: version(name) for name in
                     ("numpy", "scipy", "anndata", "scikit-learn", "joblib", "h5py")},
    }
    summaries = []
    for row in registered:
        dataset = row["evaluation_dataset"]
        if dataset not in selected:
            continue
        relative = Path(row["result_path"].split(":", 1)[1])
        path = args.artifact_root / relative
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError("input locator must stay below --artifact-root")
        if path.stat().st_size != int(row["file_size_bytes"]) or sha256(path) != row["sha256"]:
            raise ValueError("registered source hash or size differs for " + dataset)
        out = args.output_root / dataset
        if out.exists():
            completed = json.loads((out / "metadata.json").read_text())
            if completed["protocol"] != protocol or completed["source_sha256"] != row["sha256"]:
                raise ValueError("existing reproduction has a different protocol for " + dataset)
            for name, digest in completed["outputs_sha256"].items():
                if sha256(out / name) != digest:
                    raise ValueError("existing recovered output changed for " + dataset)
            if not completed["summary"]["all_outputs_reproduced"]:
                raise ValueError("existing reproduction failed its output checks for " + dataset)
            summaries.append(completed["summary"])
            print(dataset + ": resumed verified reproduction", flush=True)
            continue

        with h5py.File(path, "r") as handle:
            obs = read_elem(handle["obs"])
            var = read_elem(handle["var"])
            moments = read_elem(handle["layers/Ms"])
            velocity = read_elem(handle["layers/velocity"])
            neighbors = read_elem(handle["uns/neighbors/indices"])
            saved_phi = read_elem(handle["obsp/gv"]) if "gv" in handle["obsp"] else None
        data = ad.AnnData(X=None, shape=moments.shape, obs=obs, var=var)
        data.layers["Ms"] = moments
        data.layers["velocity"] = velocity
        data.uns["neighbors"] = {"indices": neighbors}
        print(dataset + ": initializing original GraphVelo representation", flush=True)
        selected_genes = None
        if args.mack_top:
            if "mack_score" not in data.var:
                raise ValueError("MacK scores are absent for " + dataset)
            selected_genes = data.var["mack_score"].sort_values(ascending=False).iloc[:args.mack_top].index
        model = GraphVelo(data, xkey="Ms", vkey="velocity", gene_subset=selected_genes)
        selected_mask = (np.asarray(data.var_names.isin(selected_genes)) if selected_genes is not None
                         else np.ones(data.n_vars, dtype=bool))
        selected_velocity = velocity[:, selected_mask]
        valid_mask = ~np.isnan(np.asarray(selected_velocity.sum(axis=0)).ravel())
        effective_genes = np.asarray(data.var_names, dtype=str)[selected_mask][valid_mask]

        if saved_phi is not None and not args.mack_top:
            model.T = sparse.csr_matrix(saved_phi)
            phi_source = "stored_obsp_gv"
            print(dataset + ": using stored native Phi", flush=True)
        else:
            # GraphVelo 0.1.11 train() body with explicit worker count.
            print(dataset + ": fitting original Phi with " + str(args.n_jobs) + " workers", flush=True)
            kernel = corr_kernel(model.X, model.V, model.nbrs_idx,
                                 corr_func=cos_corr, softmax_adjusted=False)
            prior = density_corrected_transition_matrix(kernel).toarray()
            del kernel
            coefficients = tangent_space_projection(
                model.X, model.V, prior, model.nbrs_idx,
                a=1, b=10, r=1, loss_func="linear", n_jobs=args.n_jobs,
            )
            del prior
            model.T = sparse.csr_matrix(coefficients)
            del coefficients
            phi_source = "refitted_native"
        if not np.isfinite(model.T.data).all():
            raise ValueError("nonfinite native coefficients for " + dataset)
        audits = {}
        reconstructed = model.project_velocity(moments)
        with h5py.File(path, "r") as handle:
            audits["velocity_gv"] = compare_arrays(
                reconstructed, read_elem(handle["layers/velocity_gv"])
            )
            del reconstructed
            for native_key in ("gv_umap", "gv_pca"):
                if native_key not in handle["obsm"]:
                    continue
                coordinate_key = "X_" + native_key[3:]
                coordinates = read_elem(handle["obsm/" + coordinate_key])
                projected = model.project_velocity(coordinates)
                audits[native_key] = compare_arrays(
                    projected, read_elem(handle["obsm/" + native_key])
                )
        native = model.T.copy()
        native.eliminate_zeros()
        normalized = native.copy()
        normalized.setdiag(0.0)
        normalized.eliminate_zeros()
        mass = np.asarray(np.abs(normalized).sum(axis=1)).ravel()
        zero = mass <= ZERO_TOLERANCE
        divisor = np.where(zero, 1.0, mass)
        normalized = normalized.multiply((1.0 / divisor)[:, None]).tocsr()
        if np.any(zero):
            normalized = normalized.tolil()
            for cell in np.flatnonzero(zero):
                normalized.rows[cell] = [int(cell)]
                normalized.data[cell] = [1.0]
            normalized = normalized.tocsr()
        if not np.allclose(np.asarray(np.abs(normalized).sum(axis=1)).ravel(), 1.0):
            raise ValueError("signed-L1 normalization failed for " + dataset)
        out.mkdir()
        sparse.save_npz(out / "phi.npz", native)
        sparse.save_npz(out / "signed_transition.npz", normalized)
        np.save(out / "cell_ids.npy", np.asarray(data.obs_names, dtype=str), allow_pickle=False)
        np.save(out / "native_directional_mass.npy", mass, allow_pickle=False)
        summary = {
            "dataset": dataset, "n_cells": data.n_obs, "phi_source": phi_source,
            "native_phi_nnz": native.nnz, "zero_directional_rows": int(zero.sum()),
            "gene_velocity_max_abs_error": audits["velocity_gv"]["max_abs_error"],
            "gene_velocity_relative_l2_error": audits["velocity_gv"]["relative_l2_error"],
            "gene_velocity_allclose": audits["velocity_gv"]["allclose"],
            "embedded_velocity_max_abs_error": max(
                v["max_abs_error"] for k, v in audits.items() if k != "velocity_gv"
            ),
            "all_outputs_reproduced": all(v["allclose"] for v in audits.values()),
            "all_outputs_entrywise_equal": all(v["entrywise_equal"] for v in audits.values()),
        }
        completed = {
            "protocol": protocol, "dataset": dataset,
            "source_locator": row["result_path"], "source_sha256": row["sha256"],
            "summary": summary, "output_comparisons": audits,
            "selected_training_genes": effective_genes.tolist(),
            "outputs_sha256": {p.name: sha256(p) for p in sorted(out.glob("*.np*"))},
        }
        (out / "metadata.json").write_text(json.dumps(completed, indent=2, sort_keys=True) + "\n")
        summaries.append(summary)
        with (args.output_root / "summary.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(summary), lineterminator="\n")
            writer.writeheader()
            writer.writerows(summaries)
        print(dataset + ": " + json.dumps(summary), flush=True)
        if not summary["all_outputs_reproduced"]:
            raise RuntimeError("native GraphVelo replay differs from saved outputs for " + dataset)
        del model, data, moments, velocity, native, normalized
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
