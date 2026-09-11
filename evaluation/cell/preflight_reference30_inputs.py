#!/usr/bin/env python3
"""Audit cell and gene compatibility before reference-30NN metric scoring."""

from __future__ import print_function

import argparse
import csv
import json
from pathlib import Path

import h5py
import numpy as np


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


def load_reference(path):
    with np.load(path, allow_pickle=False) as data:
        return {
            "obs_names": np.asarray(data["obs_names"]).astype(str),
            "var_names": np.asarray(data["var_names"]).astype(str),
            "Ms_shape": list(data["Ms"].shape),
            "Mu_shape": list(data["Mu"].shape),
            "metadata": json.loads(str(data["metadata"])),
        }


def casefold_index(names):
    index = {}
    for position, name in enumerate(names):
        index.setdefault(str(name).casefold(), []).append(position)
    return index


def gene_mapping_summary(reference_names, method_names):
    reference_index = casefold_index(reference_names)
    method_index = casefold_index(method_names)
    shared = []
    ambiguous = []
    for token, reference_positions in reference_index.items():
        method_positions = method_index.get(token, [])
        if len(reference_positions) == 1 and len(method_positions) == 1:
            shared.append(token)
        elif method_positions:
            ambiguous.append(token)
    return {
        "shared_tokens": sorted(shared),
        "shared_gene_count": len(shared),
        "ambiguous_token_count": len(ambiguous),
        "reference_duplicate_token_count": sum(
            len(positions) > 1 for positions in reference_index.values()
        ),
        "method_duplicate_token_count": sum(
            len(positions) > 1 for positions in method_index.values()
        ),
    }


def h5_layer_shape(path, layer_key):
    """Read a layer's shape from h5ad metadata without materializing its values."""
    with h5py.File(path, "r") as handle:
        layer_path = "layers/{}".format(layer_key)
        if layer_path not in handle:
            raise KeyError("missing {}".format(layer_path))
        layer = handle[layer_path]
        if isinstance(layer, h5py.Dataset):
            return list(layer.shape)
        if "shape" in layer.attrs:
            return [int(value) for value in layer.attrs["shape"]]
        raise ValueError("{} has no readable shape metadata".format(layer_path))


def h5_text_values(node):
    if isinstance(node, h5py.Dataset):
        values = node[()]
        if values.dtype.kind == "S":
            return np.asarray([value.decode("utf-8") for value in values])
        return np.asarray(values).astype(str)
    if isinstance(node, h5py.Group) and "codes" in node and "categories" in node:
        categories = h5_text_values(node["categories"])
        return categories[np.asarray(node["codes"], dtype=np.int64)]
    raise ValueError("unsupported h5ad index encoding")


def h5_axis_names(path, axis):
    with h5py.File(path, "r") as handle:
        group = handle[axis]
        index_name = group.attrs.get("_index", "_index")
        if isinstance(index_name, bytes):
            index_name = index_name.decode("utf-8")
        if index_name not in group:
            raise KeyError("{}/{} is absent".format(axis, index_name))
        return h5_text_values(group[index_name])


def inspect_row(row, reference, artifact_root):
    source_path, resolved_logical_path = resolve_source(
        artifact_root, row["result_path"]
    )
    cell_names = h5_axis_names(source_path, "obs")
    gene_names = h5_axis_names(source_path, "var")
    layer_key = row["velocity_layer_key"]
    layer_shape = h5_layer_shape(source_path, layer_key)
    reference_cells = reference["obs_names"]
    reference_cell_set = set(reference_cells)
    method_cell_set = set(cell_names)
    common_cells = reference_cell_set.intersection(method_cell_set)
    genes = gene_mapping_summary(reference["var_names"], gene_names)
    return {
            "method": row["operational_id"],
            "dataset": row["evaluation_dataset"],
            "result_path_manifest": row["result_path"],
            "result_path_resolved_logical": resolved_logical_path,
            "velocity_layer_key": layer_key,
            "n_reference_cells": int(len(reference_cells)),
            "n_method_cells": int(len(cell_names)),
            "n_common_cells": int(len(common_cells)),
            "reference_cells_missing_from_method": int(
                len(reference_cell_set.difference(method_cell_set))
            ),
            "nonreference_method_cells": int(
                len(method_cell_set.difference(reference_cell_set))
            ),
            "duplicate_method_cell_ids": int(len(cell_names) - len(method_cell_set)),
            "n_method_genes": int(len(gene_names)),
            "velocity_layer_shape": layer_shape,
            "velocity_shape_matches_adata": layer_shape == [len(cell_names), len(gene_names)],
            **genes,
            "alignment_mode_manifest": row.get("alignment_mode", ""),
            "n_common_manifest": int(row["n_common"]),
            "status": "ok",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-root", type=Path, required=True)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--reference-moments-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    rows = read_tsv(args.inputs)
    references = {}
    for dataset in sorted({row["evaluation_dataset"] for row in rows}):
        references[dataset] = load_reference(
            args.reference_moments_root / (dataset + ".npz")
        )

    records = []
    failures = []
    for position, row in enumerate(rows, start=1):
        try:
            record = inspect_row(
                row, references[row["evaluation_dataset"]], args.artifact_root
            )
        except Exception as error:
            record = {
                "method": row["operational_id"],
                "dataset": row["evaluation_dataset"],
                "status": "error",
                "error": "{}: {}".format(type(error).__name__, error),
            }
            failures.append(record)
        records.append(record)
        print("[{}/{}] {}/{}: {}".format(
            position, len(rows), record["dataset"], record["method"], record["status"]
        ), flush=True)

    per_dataset = {}
    for dataset, reference in references.items():
        successful = [
            row for row in records
            if row["dataset"] == dataset and row["status"] == "ok"
        ]
        common_tokens = set(reference["var_names"].astype(str).tolist())
        token_sets = [set(row["shared_tokens"]) for row in successful]
        all_method_tokens = set.intersection(*token_sets) if token_sets else set()
        per_dataset[dataset] = {
            "reference_moments": {
                "n_cells": int(len(reference["obs_names"])),
                "n_genes": int(len(reference["var_names"])),
                "Ms_shape": reference["Ms_shape"],
                "Mu_shape": reference["Mu_shape"],
                "expression_mode": reference["metadata"]["expression_mode"],
            },
            "n_successful_methods": len(successful),
            "genes_shared_by_all_successful_methods": len(all_method_tokens),
            "all_method_shared_gene_tokens": sorted(all_method_tokens),
            "method_cell_deviations": [
                {
                    key: row[key] for key in (
                        "method", "n_reference_cells", "n_method_cells", "n_common_cells",
                        "reference_cells_missing_from_method", "nonreference_method_cells",
                        "alignment_mode_manifest",
                    )
                }
                for row in successful
                if row["n_method_cells"] != row["n_reference_cells"]
                or row["n_common_cells"] != row["n_reference_cells"]
            ],
        }

    payload = {
        "rows": records,
        "datasets": per_dataset,
        "failures": failures,
        "n_rows": len(records),
        "n_failures": len(failures),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("wrote {} (failures={})".format(args.output, len(failures)), flush=True)
    if failures:
        raise RuntimeError("{} preflight rows failed".format(len(failures)))


if __name__ == "__main__":
    main()
