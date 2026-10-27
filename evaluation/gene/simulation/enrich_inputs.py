"""Prepare enriched copies of frozen simulation inputs without changing observations."""
import argparse
import hashlib
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def digest(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            result.update(chunk)
    return result.hexdigest()


def assert_same(left, right, label):
    if isinstance(left, dict):
        assert set(left) == set(right), label
        for key in left:
            assert_same(left[key], right[key], label + "/" + str(key))
    elif isinstance(left, pd.DataFrame):
        pd.testing.assert_frame_equal(left, right, check_exact=True)
    elif left is None:
        assert right is None, label
    else:
        a = left.toarray() if sparse.issparse(left) else np.asarray(left)
        b = right.toarray() if sparse.issparse(right) else np.asarray(right)
        assert a.shape == b.shape and a.dtype == b.dtype, label
        if a.dtype.kind in "fc":
            assert np.array_equal(a, b, equal_nan=True), label
        else:
            assert np.array_equal(a, b), label


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-root", type=Path, required=True)
    parser.add_argument("--reference-root", type=Path, required=True)
    parser.add_argument("--selected-reference-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    args.output_root.mkdir(parents=True, exist_ok=False)
    manifest_path = args.reference_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    results = []
    changed_layers = {"clean_spliced", "clean_unspliced", "velocity_s_true", "velocity_u_true"}
    for record in manifest["scenario_records"]:
        scenario = record["scenario"]
        source = args.artifact_root / record["source"]
        source_sha = digest(source)
        assert source_sha == record["source_sha256"], scenario + " input identity"
        reference_path = args.reference_root / record["output"]
        assert digest(reference_path) == record["output_sha256"], scenario + " reference identity"
        with np.load(reference_path, allow_pickle=False) as reference:
            arrays = {key: reference[key] for key in reference.files}
        original = ad.read_h5ad(source)
        assert np.array_equal(original.obs_names.astype(str), arrays["cell_ids"]), scenario + " cells"
        assert np.array_equal(original.var_names.astype(str), arrays["gene_names"]), scenario + " genes"
        assert tuple(original.shape) == tuple(arrays["clean_s"].shape), scenario + " full scope"
        for key in ["clean_s", "clean_u", "velocity_s_true", "velocity_u_true"]:
            assert arrays[key].shape == original.shape and np.isfinite(arrays[key]).all(), key
        alpha = np.asarray(original.layers["alpha_true"])
        beta = np.asarray(original.layers["beta_true"])
        gamma = np.asarray(original.layers["gamma_true"])
        np.testing.assert_array_equal(arrays["velocity_s_true"], beta * arrays["clean_u"] - gamma * arrays["clean_s"])
        np.testing.assert_array_equal(arrays["velocity_u_true"], alpha - beta * arrays["clean_u"])
        with np.load(args.selected_reference_root / record["output"], allow_pickle=False) as selected:
            np.testing.assert_array_equal(selected["cell_ids"], arrays["cell_ids"])
            columns = original.var_names.get_indexer(selected["gene_names"].astype(str))
            assert np.all(columns >= 0), scenario + " selected genes"
            for key in ["clean_s", "clean_u", "velocity_s_true", "velocity_u_true"]:
                np.testing.assert_array_equal(arrays[key][:, columns], selected[key])
            selected_genes = len(columns)
        enriched = original.copy()
        enriched.layers["clean_spliced"] = arrays["clean_s"]
        enriched.layers["clean_unspliced"] = arrays["clean_u"]
        enriched.layers["velocity_s_true"] = arrays["velocity_s_true"]
        enriched.layers["velocity_u_true"] = arrays["velocity_u_true"]
        assert "velocity_reference" not in original.uns, scenario + " existing reference metadata"
        enriched.uns["velocity_reference"] = {
            "definition": "ODE velocity at paired noise-free simulated states",
            "state_layers": ["clean_spliced", "clean_unspliced"],
            "velocity_layers": ["velocity_s_true", "velocity_u_true"],
            "coordinate_order": ["spliced", "unspliced"],
            "cellDancer_version": "1.1.7",
            "historical_input_sha256": source_sha,
            "reference_sha256": record["output_sha256"],
        }
        destination = args.output_root / record["source"]
        destination.parent.mkdir(parents=True, exist_ok=True)
        enriched.write_h5ad(destination, compression="gzip")
        saved = ad.read_h5ad(destination)
        assert_same(original.X, saved.X, "X")
        assert_same(original.obs, saved.obs, "obs")
        assert_same(original.var, saved.var, "var")
        for slot in ["obsm", "varm", "obsp", "varp"]:
            assert_same(dict(getattr(original, slot)), dict(getattr(saved, slot)), slot)
        assert set(saved.layers) == set(original.layers) | {"clean_spliced", "clean_unspliced"}
        for key in original.layers:
            if key not in changed_layers:
                assert_same(original.layers[key], saved.layers[key], "layers/" + key)
        assert set(saved.uns) == set(original.uns) | {"velocity_reference"}
        for key in original.uns:
            assert_same(original.uns[key], saved.uns[key], "uns/" + key)
        if original.raw is None:
            assert saved.raw is None
        else:
            assert_same(original.raw.X, saved.raw.X, "raw/X")
            assert_same(original.raw.var, saved.raw.var, "raw/var")
            assert_same(dict(original.raw.varm), dict(saved.raw.varm), "raw/varm")
        for layer, key in [("clean_spliced", "clean_s"), ("clean_unspliced", "clean_u"),
                           ("velocity_s_true", "velocity_s_true"), ("velocity_u_true", "velocity_u_true")]:
            np.testing.assert_array_equal(saved.layers[layer], arrays[key])
        assert digest(source) == source_sha, scenario + " original source retained"
        results.append({
            "scenario": scenario, "path": record["source"],
            "shape": list(saved.shape), "bytes": destination.stat().st_size,
            "sha256": digest(destination), "historical_input_sha256": source_sha,
            "reference_sha256": record["output_sha256"],
            "protected_existing_layers": sorted(set(original.layers) - changed_layers),
            "observations_parameters_annotations_and_order_exact": True,
            "complete_clean_reference_and_ODE_exact": True,
            "current_selected_reference_genes_exact": selected_genes,
        })
        print(scenario, saved.shape, results[-1]["sha256"], flush=True)
    report = {"schema_version": 1, "status": "validated clean-state simulation inputs",
              "reference_manifest_sha256": digest(manifest_path),
              "preparation_script_sha256": digest(__file__), "files": results}
    (args.output_root / "manifest.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
