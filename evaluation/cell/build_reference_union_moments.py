#!/usr/bin/env python3
"""Build reference Ms/Mu for the union of reference HVGs and mapped method genes.

Gene-space methods share the fixed 30-NN smoothing graph. veloAE is evaluated
in its native 100-D latent representation.
"""

from __future__ import print_function

import argparse
import csv
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np


HERE = Path(__file__).resolve().parent
DATASET_SOURCES = {
    "bonemarrow": "bonemarrow/bonemarrow.h5ad",
    "capofibro": "capo_pre/capo_pre_counts.h5ad",
    "DentateGyrus_lamanno": "DentateGyrus_lamanno/DentateGyrus_lamanno.h5ad",
    "dentategyrus": "dentategyrus/dentategyrus.h5ad",
    "forebrain": "ForebrainGlut/ForebrainGlut.loom",
    "gastrulation": "Gastrulation/Gastrulation.h5ad",
    "hematopoiesis": "hematopoiesis_dynamo/hematopoiesis_dynamo.h5ad",
    "mousebrain": "mousebrain/mousebrain_counts.h5ad",
    "neuron": "neuron/neuron_counts.h5ad",
    "organoid": "organoid/organoid_pre.h5ad",
    "pancreas": "Pancreas/Pancreas.h5ad",
    "retina": "retina/retina_full.h5ad",
    "rpe1": "rpe1/rpe1.h5ad",
    "rpe1_gen": "rpe1_gen/rpe1_gen_counts.h5ad",
}
EXPECTED_CELLS = {
    "bonemarrow": 5780,
    "capofibro": 1222,
    "DentateGyrus_lamanno": 18213,
    "dentategyrus": 2930,
    "forebrain": 1720,
    "gastrulation": 9815,
    "hematopoiesis": 1947,
    "mousebrain": 3365,
    "neuron": 3066,
    "organoid": 3831,
    "pancreas": 3696,
    "retina": 2726,
    "rpe1": 2793,
    "rpe1_gen": 3354,
}


def _read_source(path: Path, scanpy):
    if path.suffix.lower() == ".loom":
        return scanpy.read_loom(path, sparse=True)
    return scanpy.read_h5ad(path)


from reference_neighbors import REFERENCE_PARAMS, apply_reference_neighbors  # noqa: E402


HISTORICAL_MOMENT_EXPRESSION_MODE = {"rpe1": "native_X"}
GENE_ID_FIELDS = (
    "var_names",
    "query",
    "scopes",
    "gene_name_mapping",
    "symbol",
    "_id",
    "gene_ids",
    "gene_id",
    "Gene",
    "ensembl_id",
    "gene_symbol",
)


def read_tsv(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def resolve_source(artifact_root, result_path):
    prefix, separator, relative = result_path.partition(":")
    if prefix != "velocity" or not separator:
        raise ValueError("unsupported source result path: {}".format(result_path))
    relative_path = Path(relative)
    if relative_path.is_absolute() or ".." in relative_path.parts:
        raise ValueError("unsafe source result path: {}".format(result_path))
    candidate = artifact_root / "velocity" / relative_path
    if not candidate.is_file():
        raise FileNotFoundError(str(candidate))
    return candidate, result_path


def h5_strings(node):
    if isinstance(node, h5py.Dataset):
        values = node[()]
        if values.dtype.kind == "S":
            return np.asarray([value.decode("utf-8") for value in values])
        return np.asarray(values).astype(str)
    if isinstance(node, h5py.Group) and "codes" in node and "categories" in node:
        categories = h5_strings(node["categories"])
        return categories[np.asarray(node["codes"], dtype=np.int64)]
    raise ValueError("unsupported h5ad string encoding")


def method_gene_fields(path):
    with h5py.File(path, "r") as handle:
        group = handle["var"]
        index_name = group.attrs.get("_index", "_index")
        if isinstance(index_name, bytes):
            index_name = index_name.decode("utf-8")
        fields = {"var_names": h5_strings(group[index_name])}
        for field in GENE_ID_FIELDS[1:]:
            if field in group:
                try:
                    values = h5_strings(group[field])
                except ValueError:
                    continue
                if len(values) == len(fields["var_names"]):
                    fields[field] = values
    return fields


def unique_lookup(names, label):
    lookup = {}
    duplicates = set()
    for position, name in enumerate(names):
        token = str(name).casefold()
        if token in lookup:
            duplicates.add(token)
        else:
            lookup[token] = position
    if duplicates:
        raise ValueError(
            "{} has duplicate case-insensitive identifiers (for example {})".format(
                label, sorted(duplicates)[0]
            )
        )
    return lookup


def map_method_to_reference(reference_names, fields):
    reference_lookup = unique_lookup(reference_names, "reference genes")
    matched = {}
    field_counts = {}
    used_method_positions = set()
    field_order = ["var_names"] + [field for field in fields if field != "var_names"]
    for field in field_order:
        if field not in fields:
            continue
        count = 0
        for method_position, value in enumerate(fields[field]):
            token = str(value).casefold()
            reference_position = reference_lookup.get(token)
            if (
                reference_position is not None
                and reference_position not in matched
                and method_position not in used_method_positions
            ):
                matched[reference_position] = method_position
                used_method_positions.add(method_position)
                count += 1
        field_counts[field] = count
    return matched, field_counts


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_hvg_reference(dataset, root):
    with np.load(root / (dataset + ".npz"), allow_pickle=False) as data:
        return {
            "obs_names": np.asarray(data["obs_names"]).astype(str),
            "var_names": np.asarray(data["var_names"]).astype(str),
            "Ms": np.asarray(data["Ms"], dtype=np.float32),
            "Mu": np.asarray(data["Mu"], dtype=np.float32),
        }


def build_one(dataset, rows, artifact_root, neighbor_root, hvg_root, output_root, scanpy, scvelo):
    source_relative = Path("data") / DATASET_SOURCES[dataset]
    source = artifact_root / "velocity" / source_relative
    adata = _read_source(source, scanpy)
    if adata.n_obs != EXPECTED_CELLS[dataset]:
        raise ValueError("{}: unexpected source cell count".format(dataset))
    if {"spliced", "unspliced"}.difference(adata.layers):
        raise KeyError("{}: source lacks spliced/unspliced".format(dataset))
    expression_mode = HISTORICAL_MOMENT_EXPRESSION_MODE.get(dataset, "native_X")
    if expression_mode == "spliced_plus_unspliced":
        adata.X = adata.layers["spliced"] + adata.layers["unspliced"]
    elif adata.X is None:
        adata.X = adata.layers["spliced"].copy()

    # n_top_genes=None keeps every gene passing the historical count filter.
    scvelo.pp.filter_and_normalize(
        adata, min_shared_counts=20, n_top_genes=None, log=True
    )
    reference_names = np.asarray(adata.var_names).astype(str)
    reference_lookup = unique_lookup(reference_names, "filtered reference genes")
    hvg = load_hvg_reference(dataset, hvg_root)
    hvg_positions = [reference_lookup[str(name).casefold()] for name in hvg["var_names"]]

    union_positions = set(hvg_positions)
    method_coverage = {}
    unmapped_gene_space_methods = []
    rows_by_method = {row["operational_id"]: row for row in rows}
    for row in rows:
        method = row["operational_id"]
        result_path, logical_path = resolve_source(artifact_root, row["result_path"])
        fields = method_gene_fields(result_path)
        if method == "graphvelo_dynamo":
            parent_row = rows_by_method["dynamo"]
            parent_path, _ = resolve_source(artifact_root, parent_row["result_path"])
            parent_fields = method_gene_fields(parent_path)
            parent_lookup = unique_lookup(
                parent_fields["var_names"], "Dynamo parent var_names"
            )
            parent_positions = np.asarray(
                [parent_lookup.get(str(name).casefold(), -1) for name in fields["var_names"]],
                dtype=np.int64,
            )
            for field, values in parent_fields.items():
                if field != "var_names":
                    inherited = np.full(len(parent_positions), "", dtype=object)
                    valid = parent_positions >= 0
                    inherited[valid] = values[parent_positions[valid]]
                    fields["parent_dynamo_" + field] = inherited
        mapped, field_counts = map_method_to_reference(reference_names, fields)
        is_latent_veloae = method == "veloae" and not mapped
        if not is_latent_veloae:
            union_positions.update(mapped)
            if not mapped:
                unmapped_gene_space_methods.append(method)
        method_coverage[method] = {
            "result_path": logical_path,
            "mapped_reference_gene_count": len(mapped),
            "mapping_field_counts": field_counts,
            "coordinate_space": "method_native_latent" if is_latent_veloae else "reference_gene_space",
        }
    if unmapped_gene_space_methods:
        raise RuntimeError(
            "{}: non-latent methods without a reference-gene mapping: {}".format(
                dataset, ", ".join(unmapped_gene_space_methods)
            )
        )

    selected_positions = np.asarray(sorted(union_positions), dtype=np.int64)
    adata = adata[:, selected_positions].copy()
    applied, detail = apply_reference_neighbors(
        adata, dataset, output_root=neighbor_root
    )
    if not applied:
        raise RuntimeError("{}: {}".format(dataset, detail))
    scvelo.pp.moments(adata, n_neighbors=REFERENCE_PARAMS["n_neighbors"])

    selected_lookup = unique_lookup(np.asarray(adata.var_names).astype(str), "union genes")
    hvg_union_positions = np.asarray(
        [selected_lookup[str(name).casefold()] for name in hvg["var_names"]], dtype=np.int64
    )
    validation = {
        "obs_names_exact": bool(np.array_equal(np.asarray(adata.obs_names).astype(str), hvg["obs_names"])),
        "hvg_names_exact": bool(np.array_equal(np.asarray(adata.var_names).astype(str)[hvg_union_positions], hvg["var_names"])),
        "Ms_hvg_exact": bool(np.array_equal(adata.layers["Ms"][:, hvg_union_positions], hvg["Ms"])),
        "Mu_hvg_exact": bool(np.array_equal(adata.layers["Mu"][:, hvg_union_positions], hvg["Mu"])),
    }
    if not all(validation.values()):
        failed = [key for key, value in validation.items() if not value]
        raise RuntimeError("{}: reference-moment validation failed: {}".format(dataset, failed))

    output_root.mkdir(parents=True, exist_ok=True)
    output = output_root / (dataset + ".npz")
    metadata = {
        "dataset": dataset,
        "source": "velocity:" + source_relative.as_posix(),
        "source_sha256": sha256(source),
        "expression_mode": expression_mode,
        "preprocessing": "scvelo_0.2.5/filter_and_normalize(min_shared_counts=20,n_top_genes=None,log=True)",
        "smoothing_graph": "reference_neighbors/scvelo_0.2.5/fixed_30nn",
        "selected_gene_rule": "historical_2000_hvgs UNION mapped_gene_space_method_outputs",
        "historical_hvg_validation": validation,
        "method_coverage": method_coverage,
    }
    np.savez_compressed(
        output,
        obs_names=np.asarray(adata.obs_names).astype(str),
        var_names=np.asarray(adata.var_names).astype(str),
        Ms=np.asarray(adata.layers["Ms"], dtype=np.float32),
        Mu=np.asarray(adata.layers["Mu"], dtype=np.float32),
        metadata=np.asarray(json.dumps(metadata, sort_keys=True)),
    )
    return {
        "dataset": dataset,
        "output": str(output),
        "n_cells": int(adata.n_obs),
        "n_union_genes": int(adata.n_vars),
        "status": "ok",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-root", type=Path, required=True)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument(
        "--reference-neighbors-root",
        type=Path,
        default=HERE / "reference_neighbors",
        help="fixed reference-30 assets; defaults to the repository copy",
    )
    parser.add_argument("--hvg-reference-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--datasets", nargs="*", choices=sorted(DATASET_SOURCES))
    args = parser.parse_args()

    import scanpy as sc
    import scvelo as scv

    rows = read_tsv(args.inputs)
    selected = args.datasets or sorted(DATASET_SOURCES)
    manifest_path = args.output_root / "manifest.json"
    if manifest_path.is_file():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        manifest["scvelo_version"] = scv.__version__
        manifest.setdefault("datasets", {})
    else:
        manifest = {"datasets": {}, "scvelo_version": scv.__version__}
    failures = []
    for dataset in selected:
        dataset_rows = [row for row in rows if row["evaluation_dataset"] == dataset]
        if len(dataset_rows) != 20:
            raise ValueError("{}: expected 20 method rows, found {}".format(dataset, len(dataset_rows)))
        try:
            result = build_one(
                dataset, dataset_rows, args.artifact_root,
                args.reference_neighbors_root, args.hvg_reference_root,
                args.output_root, sc, scv,
            )
        except Exception as error:
            result = {"dataset": dataset, "status": "error", "error": "{}: {}".format(type(error).__name__, error)}
            failures.append(dataset)
        manifest["datasets"][dataset] = result
        print(result, flush=True)
    args.output_root.mkdir(parents=True, exist_ok=True)
    (args.output_root / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    if failures:
        raise RuntimeError("failed datasets: {}".format(", ".join(failures)))


if __name__ == "__main__":
    main()
