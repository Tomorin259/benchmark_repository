#!/usr/bin/env python3
"""Recompute the selected normalized GeoMedian-5 fields, scores and S5 ablations."""
import argparse
import csv
import os
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from importlib.metadata import version
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
if __package__ in {None, ""}:
    sys.path.insert(0, str(ROOT))

from ensemble.benchmark_normalization import METHODS, aggregate_rows, normalize_offdiagonal
from ensemble.core import DEFAULT_CONFIG, _canonical_dataset, load_config, row_l1_norm, sha256, validate_scvelo_runtime
from ensemble.pareto_members import pareto_front, score_fields
from ensemble.reproject import project_graphvelo
from ensemble.supfig10_data import write_panel_data
from ensemble.pareto_bootstrap import SOURCE_NAMES as PANEL_SOURCES

METRICS = ("cbdir", "batc")
CONFIGURATIONS = {
    "full5": ("GeoMedian-5 (full)", METHODS),
    "pair2": ("TIvelo + veloVI (2 members)", METHODS[:2]),
    **{
        "without_" + method.lower().removesuffix("_scv"): (
            "w/o " + method + (" (GeoMedian-4)" if method == "GraphVelo_scv" else ""),
            tuple(other for other in METHODS if other != method),
        ) for method in METHODS
    },
    "mean5": ("Mean-5 (same 5 members)", METHODS),
}


def replay_dataset(dataset, *, input_root, projections, specs, records):
    folder = input_root / dataset
    record = records[dataset]
    n = record["n_cells"]
    names = {"geometry.npz", *("operator_input_" + m + ".npz" for m in METHODS)}
    if set(record["sha256"]) != names:
        raise ValueError(f"{dataset}: expected geometry and five member checksums")
    for name in names:
        if sha256(folder / name) != record["sha256"][name]:
            raise ValueError(f"{dataset}/{name}: checksum mismatch")
    bases = ("pca", "umap") if dataset == "forebrain" else ("umap",)
    with np.load(folder / "geometry.npz", allow_pickle=False) as archive:
        geometry = {key: archive[key] for key in ("cell_ids", "labels", "neighbors", *("X_" + b for b in bases))}
    cells, labels, neighbors = (geometry[key] for key in ("cell_ids", "labels", "neighbors"))
    if (cells.shape != (n,) or cells.dtype.kind != "U" or len(set(cells)) != n
            or labels.shape != (n,) or labels.dtype.kind != "U" or "" in cells or "" in labels):
        raise ValueError(f"{dataset}: invalid cell IDs or label axis")
    if (neighbors.shape != (n, 30) or neighbors.dtype.kind not in "iu"
            or neighbors.min() < 0 or neighbors.max() >= n
            or not np.array_equal(neighbors[:, 0], np.arange(n))):
        raise ValueError(f"{dataset}: invalid frozen neighbor indices")
    spec = specs[dataset]
    plot = ad.AnnData(obs=pd.DataFrame({spec["label"]: labels}, index=cells))
    plot.uns["neighbors"] = {"indices": neighbors}
    for basis in bases:
        coords = geometry["X_" + basis]
        if coords.shape != (n, 30 if basis == "pca" else 2) or not np.isfinite(coords).all():
            raise ValueError(f"{dataset}: invalid {basis} coordinates")
        plot.obsm["X_" + basis] = coords

    normalized = {}
    for method in METHODS:
        matrix = sparse.load_npz(folder / f"operator_input_{method}.npz").tocsr()
        matrix.check_format(full_check=True)
        if (matrix.shape != (n, n) or not np.isfinite(matrix.data).all()
                or np.any(matrix.diagonal()) or np.any(row_l1_norm(matrix) > 1 + 1e-10)):
            raise ValueError(f"{dataset}/{method}: invalid frozen signed weights")
        normalized[method] = normalize_offdiagonal(matrix)
        if method == "GraphVelo_scv" and np.max(abs((normalized[method] - matrix).data), initial=0) > 1e-12:
            raise ValueError(f"{dataset}: native GraphVelo Phi must already have unit active row mass")
    member_fields = {b: {m: project_graphvelo(w, geometry["X_" + b])
                         for m, w in normalized.items()} for b in bases}
    rows = []
    for key, (_, members) in CONFIGURATIONS.items():
        group = [normalized[m] for m in members]
        if key == "mean5":
            operator = sum(group[1:], group[0]) * (1 / len(group))
            coefficients = np.full((n, len(group)), 1 / len(group))
        else:
            operator, diagnostic, coefficients = aggregate_rows(group)
            if diagnostic["n_converged_rows"] != n:
                raise ValueError(f"{dataset}/{key}: geometric median did not converge")
        if (np.any(operator.diagonal()) or np.any(row_l1_norm(operator) > 1 + 1e-10)
                or coefficients.min() < 0 or not np.allclose(coefficients.sum(1), 1, rtol=0, atol=1e-12)
                or (key == "pair2" and not np.allclose(coefficients, .5, rtol=0, atol=1e-12))):
            raise ValueError(f"{dataset}/{key}: invalid aggregate or coefficients")
        fields = {b: project_graphvelo(operator, geometry["X_" + b]) for b in bases}
        error = max(float(np.max(abs(field - sum(coefficients[:, i, None] * member_fields[b][m]
                    for i, m in enumerate(members))))) for b, field in fields.items())
        if not np.isfinite(error) or error > 1e-9:
            raise ValueError(f"{dataset}/{key}: projection reconstruction failed")
        scored = score_fields(plot, fields, dataset, spec)
        if not all(np.isfinite(scored[m]) for m in METRICS):
            raise ValueError(f"{dataset}/{key}: nonfinite score")
        rows.append(dict(dataset=dataset, n_cells=n, id=key, members="|".join(members),
                         aggregation="mean" if key == "mean5" else "row_geometric_median",
                         cbdir=scored["cbdir"], batc=scored["batc"],
                         cbdir_evaluable_source_cells=scored["cbdir_edge_evaluable_source_cells"],
                         batc_evaluable_source_cells=scored["batc_evaluable_source_cells"],
                         projection_reconstruction_error=error))
        if key == "full5":
            np.savez_compressed(projections / (dataset + ".npz"), **geometry,
                                **{"ensemble_velocity_" + b: v for b, v in fields.items()})
    print(f"{dataset}: scored eight configurations", flush=True)
    return rows


def run(input_root, output, jobs):
    validate_scvelo_runtime()
    config = load_config()
    manifest_path = input_root / "manifest.json"
    with (ROOT / "artifacts.tsv").open(newline="") as stream:
        registered = [row["sha256"] for row in csv.DictReader(stream, delimiter="\t")
                      if row["artifact_id"] == "ensemble-geomedian5-frozen-input-manifest"]
    if registered != [sha256(manifest_path)]:
        raise ValueError("Frozen input manifest does not match the registered GeoMedian-5 release")
    manifest = json.loads(manifest_path.read_text())
    specs = config["datasets"]
    records = manifest["datasets"]
    if (manifest["schema_version"] != 1 or manifest["methods"] != list(METHODS)
            or manifest["config_sha256"] != sha256(DEFAULT_CONFIG) or set(records) != set(specs)):
        raise ValueError("Frozen-input schema, member order, configuration or dataset set differs")
    dataset_table = ROOT / "data/datasets.tsv"
    mapping_table = ROOT / "evaluation/cell/inputs.tsv"
    registry = pd.read_csv(dataset_table, sep="\t")
    aliases = pd.read_csv(mapping_table, sep="\t", usecols=["dataset_id", "evaluation_dataset"]).drop_duplicates()
    aliases["dataset"] = aliases.evaluation_dataset.map(_canonical_dataset)
    registry = registry.merge(aliases[["dataset_id", "dataset"]].drop_duplicates(), on="dataset_id",
                              validate="one_to_one", sort=False)
    registry = registry[registry.dataset.isin(specs)].set_index("dataset", verify_integrity=True)
    if set(registry.index) != set(specs) or any(records[d]["n_cells"] != registry.at[d, "n_obs"] for d in specs):
        raise ValueError("Dataset registry identity or cell count differs from frozen input")
    projections = output / "projections"
    projections.mkdir()
    compute = partial(replay_dataset, input_root=input_root, projections=projections, specs=specs, records=records)
    with ProcessPoolExecutor(max_workers=jobs) as pool:
        scores = pd.DataFrame([row for block in pool.map(compute, specs) for row in block])
    selected = scores[scores.id == "full5"].set_index("dataset")
    summary = []
    for m in METRICS:
        scores["delta_" + m] = scores[m] - scores.dataset.map(selected[m])
        summary.append(dict(metric=m, n_datasets=len(selected), n_cells=selected.n_cells.sum(),
                            mean=selected[m].mean(), median=selected[m].median()))
    scores.to_csv(output / "dataset_scores.csv", index=False)
    pd.DataFrame(summary).to_csv(output / "summary.csv", index=False)
    selected = selected.assign(variant="normalized_geomedian5", label="GeoMedian-5",
                               cbdir_basis=["pca" if d == "forebrain" else "umap" for d in selected.index],
                               batc_basis="umap")
    selected[["variant", "label", "n_cells", *METRICS, "cbdir_basis", "batc_basis"]].to_csv(output / "scores.csv")
    rows = []
    for key, (label, _) in CONFIGURATIONS.items():
        group = scores[scores.id == key]
        rec = dict(id=key, configuration=label, status="ACTUAL_RUN", n_datasets=len(group))
        for m in METRICS:
            delta = group["delta_" + m]
            rec.update({m: group[m].mean(), "delta_" + m: delta.mean(),
                        m + "_higher_than_full5": (delta > 1e-12).sum(),
                        m + "_lower_than_full5": (delta < -1e-12).sum()})
        rows.append(rec)
    macro = pd.DataFrame(rows)
    macro["pareto_among_eight"] = pareto_front(macro[list(METRICS)])[0]
    macro.to_csv(output / "Table_S5_A.csv", index=False)
    table = selected.loc[registry.index, list(METRICS)].copy()
    table.insert(0, "display_name", registry.paper_name)
    table.loc["rpe1", "display_name"] = "RPE1 (scEU-seq)"
    table.loc["rpe1_gen", "display_name"] = "RPE1 (VeloCycle)"
    table.loc["macro_average"] = ["Macro-average", *(selected[m].mean() for m in METRICS)]
    table.to_csv(output / "Table_S5_B.csv")
    write_panel_data(output / "Table_S5_B.csv", output / "s10")
    code = ["__main__.py", "supfig10_data.py", "pareto_bootstrap.py", "replay_selected.py", "benchmark_normalization.py", "core.py", "geometric_median.py",
            "reproject.py", "pareto_members.py", "metrics.py"]
    protocol = dict(selected="full5", methods=list(METHODS), member_operators=["T", "T", "T", "T", "native Phi"],
                    configurations={k: list(m) for k, (_, m) in CONFIGURATIONS.items()},
                    normalization="Off-diagonal signed L1 per member row; mass <= 1e-12 becomes zero",
                    aggregation="Equal-weight Euclidean geometric median on union support; mean5 uses arithmetic mean",
                    projection="Raw weighted displacements; no centering or post-aggregation normalization",
                    evaluation="Forebrain PCA30 CBDir and UMAP2 BATC; UMAP2 for both metrics elsewhere",
                    summary="Equal dataset weights; signed delta = variant minus full5; comparison tolerance 1e-12",
                    randomness="Deterministic S5; S10a uses 10000 paired dataset bootstraps with seed 20261109", jobs=jobs, input_manifest_sha256=sha256(manifest_path),
                    config_sha256=manifest["config_sha256"],
                    registry_sha256={str(p.relative_to(ROOT)): sha256(p) for p in (dataset_table, mapping_table)},
                    panel_source_sha256={name: sha256(ROOT / name) for name in PANEL_SOURCES},
                    threads={name: os.environ.get(name) for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS")},
                    code_sha256={"ensemble/" + f: sha256(ROOT / "ensemble" / f) for f in code},
                    versions={p: version(p) for p in ("numpy", "scipy", "pandas", "anndata", "scvelo", "scikit-learn")})
    (output / "protocol.json").write_text(json.dumps(protocol, indent=2) + "\n")
    metadata = dict(status="complete", selected="full5", protocol="protocol.json", n_cells=int(selected.n_cells.sum()),
                    n_datasets=len(selected), n_configurations=len(macro), n_dataset_configurations=len(scores),
                    outputs_sha256={str(p.relative_to(output)): sha256(p) for p in sorted(output.rglob("*")) if p.is_file()})
    (output / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
        epilog="For the historical four-member H5AD workflow use ensemble/benchmark.py --legacy-four --help.")
    parser.add_argument("--input-root", type=Path, default=ROOT / "ensemble/frozen_inputs")
    parser.add_argument("--output-root", required=True, type=Path, help="new directory for the complete replay")
    parser.add_argument("--jobs", type=int, default=1)
    args = parser.parse_args(argv)
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    output = args.output_root.resolve()
    if output.is_relative_to(ROOT / "figures/assembly") or output.is_relative_to(args.input_root.resolve()):
        parser.error("choose a new output outside the input bundle and figure assembly")
    if output.exists():
        parser.error("--output-root already exists")
    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix="." + output.name + "-", dir=output.parent) as temp:
        stage = Path(temp) / "bundle"
        stage.mkdir()
        run(args.input_root, stage, args.jobs)
        output.mkdir()  # Reserve the name exclusively before publishing by rename.
        stage.replace(output)


if __name__ == "__main__":
    main()
