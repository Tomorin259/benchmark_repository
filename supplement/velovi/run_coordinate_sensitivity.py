#!/usr/bin/env python3
"""Recompute veloVI Table S6 from hash-identified, fixed-reference inputs."""

import argparse
import csv
import hashlib
import importlib.metadata
import json
import statistics
import sys
from pathlib import Path

CONDITIONS = (
    "native_ms_native_velocity",
    "reference_ms_native_velocity",
    "reference_ms_affine_velocity",
)
TABLE_FIELDS = (
    "dataset", "A_native_Ms_V", "B_reference_Ms_V",
    "C_reference_Ms_affine_V", "C_minus_B",
)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def write_table_s6(rows, path):
    """Generate the unweighted dataset-level table from the 42 source rows."""
    by_dataset = {}
    for row in rows:
        values = by_dataset.setdefault(row["dataset"], {})
        condition = row["condition"]
        if row["status"] != "ok" or condition in values:
            raise ValueError("invalid or duplicate dataset-condition result")
        value = float(row["transcosine"])
        if not float("-inf") < value < float("inf"):
            raise ValueError("nonfinite TransCosine")
        values[condition] = value
    if len(by_dataset) != 14:
        raise ValueError("Table S6 requires exactly 14 datasets")
    output = []
    for dataset, values in by_dataset.items():
        if set(values) != set(CONDITIONS):
            raise ValueError(f"incomplete conditions for {dataset}")
        a, b, c = (values[condition] for condition in CONDITIONS)
        output.append(dict(zip(TABLE_FIELDS, (dataset, a, b, c, c - b))))
    macro = {"dataset": "Macro-average"}
    for key in TABLE_FIELDS[1:]:
        macro[key] = statistics.fmean(row[key] for row in output)
    output.append(macro)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=TABLE_FIELDS)
        writer.writeheader()
        writer.writerows(output)


def align_reference_ms(reference_ms, reference_names, scored_names):
    """Reorder reference rows into the scored-cell order, never by position."""
    import numpy as np

    if len(set(reference_names)) != len(reference_names):
        raise ValueError("duplicate reference cell identifiers")
    if len(set(scored_names)) != len(scored_names):
        raise ValueError("duplicate scored cell identifiers")
    if set(reference_names) != set(scored_names):
        raise ValueError("reference and scored cell sets differ")
    positions = {name: index for index, name in enumerate(reference_names)}
    return reference_ms[np.asarray([positions[name] for name in scored_names])]


def fit_affine(native_ms, reference_ms):
    """OLS with intercept; retain the historical convention for constant genes."""
    import numpy as np

    native_mean = native_ms.mean(axis=0)
    reference_mean = reference_ms.mean(axis=0)
    centered_native = native_ms - native_mean
    variance = np.sum(centered_native * centered_native, axis=0)
    slope = np.divide(
        np.sum(centered_native * (reference_ms - reference_mean), axis=0),
        variance, out=np.ones(native_ms.shape[1]), where=variance > 1e-15,
    )
    intercept = reference_mean - slope * native_mean
    fitted = native_ms * slope + intercept
    total = np.sum((reference_ms - reference_mean) ** 2, axis=0)
    r2 = np.ones(native_ms.shape[1])
    valid = total > 1e-15
    r2[valid] = 1 - np.sum((reference_ms - fitted)[:, valid] ** 2, axis=0) / total[valid]
    return slope, intercept, r2


def mean_metric(raw_values):
    """Mean finite source-cell scores per edge, then unweighted mean over edges."""
    import numpy as np

    values = []
    for raw in raw_values.values():
        finite = np.asarray(raw, dtype=float)
        finite = finite[np.isfinite(finite)]
        if finite.size:
            values.append(float(np.mean(finite)))
    return float(np.mean(values)) if values else float("nan")


def main():
    module = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=module / "results/manifest.json")
    parser.add_argument("--artifact-root", type=Path)
    parser.add_argument("--reference-moments-root", type=Path)
    parser.add_argument("--reference-neighbors-root", type=Path,
                        default=module.parents[1] / "evaluation/cell/reference_neighbors")
    parser.add_argument("--output-root", type=Path, default=module / "results")
    parser.add_argument("--table-only", action="store_true",
                        help="Regenerate table_s6.csv from the archived CSV without scientific dependencies")
    args = parser.parse_args()
    args.output_root.mkdir(parents=True, exist_ok=True)
    if args.table_only:
        with (module / "results/full_dataset_three_condition.csv").open(newline="") as handle:
            write_table_s6(list(csv.DictReader(handle)), args.output_root / "table_s6.csv")
        return
    if args.artifact_root is None or args.reference_moments_root is None:
        parser.error("--artifact-root and --reference-moments-root are required for recomputation")

    import anndata as ad
    import numpy as np
    import scvelo as scv
    from scipy import sparse
    from tivelo.utils import metrics
    from tivelo.utils.metrics import cross_boundary_correctness2, cross_boundary_scvelo_probs

    manifest = json.loads(args.manifest.read_text())
    software = {name: importlib.metadata.version(name)
                for name in ("scvelo", "tivelo", "numpy", "scipy", "anndata", "pandas")}
    if software["scvelo"] != "0.3.1":
        raise ValueError("Table S6 requires scvelo==0.3.1")
    metric_hash = sha256(Path(metrics.__file__))
    if metric_hash != manifest["software"]["tivelo_metrics_sha256"]:
        raise ValueError("TIvelo metric source differs from the verified implementation")
    rows = []
    for item in manifest["inputs"]:
        dataset = item["dataset"]
        paths = {
            "h5ad": args.artifact_root / item["result_path"].removeprefix("velocity:"),
            "moments": args.reference_moments_root / f"{dataset}.npz",
            "neighbors": args.reference_neighbors_root / f"{dataset}.npz",
        }
        for name, path in paths.items():
            if sha256(path) != item[f"{name}_sha256"]:
                raise ValueError(f"{dataset}: {name} SHA-256 differs from Table S6 input")
        print(f"START {dataset}", flush=True)
        adata = ad.read_h5ad(paths["h5ad"])
        if not adata.obs_names.is_unique:
            raise ValueError(f"{dataset}: duplicate H5AD cell identifiers")
        with np.load(paths["neighbors"], allow_pickle=False) as npz:
            names = np.asarray(npz["obs_names"]).astype(str)
            adata = adata[names].copy()
            for key in ("distances", "connectivities"):
                shape = tuple(int(value) for value in npz[f"{key}_shape"])
                adata.obsp[key] = sparse.csr_matrix(
                    (npz[f"{key}_data"], npz[f"{key}_indices"], npz[f"{key}_indptr"]),
                    shape=shape,
                )
            adata.uns["neighbors"] = {
                "indices": np.asarray(npz["indices"], dtype=np.int64),
                "distances_key": "distances", "connectivities_key": "connectivities",
                "params": {"n_neighbors": 30},
            }
        with np.load(paths["moments"], allow_pickle=False) as npz:
            reference_names = np.asarray(npz["obs_names"]).astype(str)
            reference_genes = np.asarray(npz["var_names"]).astype(str)
            reference_ms = np.asarray(npz["Ms"], dtype=np.float64)
        reference_ms = align_reference_ms(reference_ms, reference_names, names)
        # The historical common-gene contract uses case-insensitive identifiers.
        lookup = {name.casefold(): index for index, name in enumerate(reference_genes)}
        if len(lookup) != len(reference_genes):
            raise ValueError(f"{dataset}: ambiguous case-insensitive reference gene identifiers")
        method_names = [str(name).casefold() for name in adata.var_names]
        if len(set(method_names)) != len(method_names):
            raise ValueError(f"{dataset}: ambiguous case-insensitive method gene identifiers")
        method_genes = [index for index, name in enumerate(method_names) if name in lookup]
        if not method_genes:
            raise ValueError(f"{dataset}: no common genes")
        adata = adata[:, method_genes].copy()
        reference_ms = reference_ms[:, [lookup[method_names[index]] for index in method_genes]]
        native_ms = adata.layers["Ms"]
        native_v = adata.layers["velocity"]
        native_ms = (native_ms.toarray() if sparse.issparse(native_ms) else np.asarray(native_ms)).astype(np.float64)
        native_v = (native_v.toarray() if sparse.issparse(native_v) else np.asarray(native_v)).astype(np.float64)
        if not all(np.isfinite(x).all() for x in (native_ms, native_v, reference_ms)):
            raise ValueError(f"{dataset}: nonfinite input matrices")
        slope, intercept, r2 = fit_affine(native_ms, reference_ms)
        adata.obs["_velovi_s6_label"] = np.asarray(adata.obs[item["label"]])
        edges = [tuple(edge) for edge in item["cluster_edges"]]
        for condition, ms, velocity in zip(
            CONDITIONS, (native_ms, reference_ms, reference_ms),
            (native_v, native_v, native_v * slope[None, :]),
        ):
            adata.layers["Ms"], adata.layers["velocity"] = ms, velocity
            for key in ("velocity_graph", "velocity_graph_neg", "Velocity_graph", "Velocity_graph_neg"):
                adata.uns.pop(key, None)
            scv.tl.velocity_graph(adata, vkey="velocity", xkey="Ms", n_neighbors=30,
                                  n_recurse_neighbors=1, mode_neighbors="distances", n_jobs=8)
            rows.append({
                "dataset": dataset, "condition": condition, "status": "ok",
                "highdim_cbdir": mean_metric(cross_boundary_correctness2(
                    adata, "_velovi_s6_label", "velocity", edges, return_raw=True)),
                "transcosine": mean_metric(cross_boundary_scvelo_probs(
                    adata, "_velovi_s6_label", edges, "velocity_graph", return_raw=True)),
                "n_cells": adata.n_obs, "n_common_genes": len(method_genes),
                "n_edges": len(edges), "n_neighbors": 30,
                "affine_median_slope": float(np.median(slope)),
                "affine_mean_slope": float(np.mean(slope)),
                "affine_slope_q05": float(np.quantile(slope, .05)),
                "affine_slope_q95": float(np.quantile(slope, .95)),
                "affine_median_r2": float(np.median(r2)),
                "affine_negative_slope_fraction": float(np.mean(slope < 0)),
                "velocity_norm_ratio": float(np.mean(
                    np.linalg.norm(native_v * slope[None, :], axis=1)
                    / np.maximum(np.linalg.norm(native_v, axis=1), 1e-12))),
                **{f"{name}_sha256": item[f"{name}_sha256"] for name in paths},
            })
            print(f"DONE {dataset} {condition} {rows[-1]['transcosine']}", flush=True)

    write_table_s6(rows, args.output_root / "table_s6.csv")
    with (args.output_root / "full_dataset_three_condition.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    manifest["result_origin"] = "recomputed with public diagnostic script"
    manifest["software"] = {**software, "python": sys.version.split()[0],
                            "tivelo_metrics_sha256": metric_hash}
    manifest["script_sha256"] = sha256(Path(__file__).resolve())
    manifest["outputs"] = {
        name: sha256(args.output_root / name)
        for name in ("full_dataset_three_condition.csv", "table_s6.csv")
    }
    (args.output_root / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
