"""Publish validated compact MultiVelo diagnostics from identified raw runs.

Run on the execution host. Large native H5ADs remain external; only exact
plot inputs, compact results and portable provenance are written here.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import shutil
import sys


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle))


def write(path, rows):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def publish_materials(source_root, output):
    materials = json.loads((source_root / "manifest.json").read_text())
    assert materials["script_sha256"] == sha(output / "plot_materials.py")
    assert materials["input_hashes_unchanged"]
    assert all(sha(output / name) == expected for name, expected in materials["inputs"].items())
    published = {}
    for name, expected in materials["outputs"].items():
        assert Path(name).name == name
        source = source_root / name
        assert sha(source) == expected
        folder = "tables" if name.startswith("controls_table_preview.") or name == "tableS9.tex" else "figures"
        target = output / folder / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        published[str(target.relative_to(output))] = sha(target)
    materials["outputs"] = published
    materials["scope"] = "Source materials for manuscript Supplementary Figure S15 A/B and Table S9"
    (output / "materials_manifest.json").write_text(json.dumps(materials, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--diagnostics-root", type=Path)
    parser.add_argument("--artifact-root", type=Path)
    parser.add_argument("--output-root", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--materials-dir", type=Path, help="Accepted output of plot_materials.py")
    parser.add_argument("--materials-only", action="store_true", help="Publish verified rendering without regenerating frozen scientific inputs or results")
    args = parser.parse_args()
    output = args.output_root.resolve()
    if args.materials_only:
        if args.materials_dir is None:
            parser.error("--materials-only requires --materials-dir")
        publish_materials(args.materials_dir, output)
        return
    if args.diagnostics_root is None or args.artifact_root is None:
        parser.error("Full result import requires --diagnostics-root and --artifact-root")
    root = args.diagnostics_root.resolve()
    repo = Path(__file__).resolve().parents[2]
    for name in ["results", "inputs"]:
        (output / name).mkdir(parents=True, exist_ok=True)

    def portable(value):
        if isinstance(value, dict):
            return {portable(k): portable(v) for k, v in value.items()}
        if isinstance(value, list):
            return [portable(x) for x in value]
        if isinstance(value, str) and value.startswith("/"):
            path = Path(value)
            for base, prefix in [(repo, "repository:"), (args.artifact_root.resolve(), "artifact-root:"), (Path(sys.prefix), "environment:")]:
                try:
                    return prefix + path.relative_to(base).as_posix()
                except ValueError:
                    pass
            return "external:" + path.name
        return value

    copies = {
        "scvelo_controls.csv": root / "scvelo_saved_fit_controls_run5/gene_controls.csv",
        "scvelo_reproduction_checks.csv": root / "scvelo_saved_fit_controls_run5/reproduction_checks.csv",
        "phase_rule_inputs.csv": root / "phase_rule_inputs.csv",
        "candidate_gate_branch_scores.csv": root / "candidate_gate_branching/branch_scores.csv",
        "candidate_gate_baseline_reproduction.csv": root / "candidate_gate_branching/baseline_reproduction.csv",
        "candidate_gate_pilot_cell_states.csv": root / "candidate_gate_pilot/cell_velocity_states.csv",
        "constant_direction_baselines.csv": root / "constant_direction_baselines.csv",
    }
    for name, source in copies.items():
        shutil.copyfile(source, output / "results" / name)
    for name in ["initialization_audit.json", "patch_validation.json", "candidate_gate_only.patch"]:
        shutil.copyfile(root / "candidate_gate_branching" / name, output / "results" / name)
    write(output / "results/input_identity_checks.csv", [r for r in read(root / "frozen_output_ablation/input_and_reconstruction_checks.csv") if r["layer"] != "MultiVelo_full_reconstruction"])
    phase_rows = read(root / "assignment_preview_v3/cell_phase_direction_audit.csv")
    write(output / "results/scvelo_phase_assignment_case.csv", [dict(cell_id=r["cell_id"], gene="BoostALL_5", bin=r["branch"] + "_" + r["progression_bin"], final_state=r["scvelo_final_state"]) for r in phase_rows])

    refits = read(root / "all_gene_paired_phase_refits/paired_refit_scores.csv")
    gate = read(root / "candidate_gate_branching/paired_refit_scores.csv")
    smooth = read(root / "velocity_smoothing_ablation/paired_refit_scores.csv")
    defaults = {(r["dataset"], r["gene"]): r for r in refits if r["condition"] == "refit_default"}
    fields = ["cbdir", "cossim", "alpha", "beta", "gamma", "switch", "unspliced_scaling", "fit_loss"]
    for row in gate:
        if row["condition"] == "refit_default":
            assert all(float(row[k]) == float(defaults[row["dataset"], row["gene"]][k]) for k in fields)
    common = list(refits[0])
    common.remove("wall_seconds")
    combined = [{key: row[key] for key in common} for row in refits + [r for r in gate if r["condition"] == "on_allow_repression"] + smooth]
    assert len(combined) == 250
    assert len({(r["dataset"], r["gene"], r["condition"]) for r in combined}) == 250
    write(output / "results/multivelo_refits.csv", combined)

    import anndata as ad
    import numpy as np
    original = ad.read_h5ad(root / "candidate_gate_branching/branching_refit_default.h5ad")[:, ["BoostALL_5"]].copy()
    opened = ad.read_h5ad(root / "candidate_gate_branching/branching_on_allow_repression.h5ad")[original.obs_names, original.var_names].copy()
    for key in ["Ms", "Mu", "velocity_s_true", "velocity_u_true"]:
        np.testing.assert_array_equal(original.layers[key], opened.layers[key])
    labels = original.obs["bin"].astype(str).to_numpy()
    arrays = {"cell_ids": original.obs_names.to_numpy(dtype=str), "bins": np.array([int(x.rsplit("_", 1)[1]) for x in labels]), "branches": np.array([x.rsplit("_", 1)[0] for x in labels])}
    for key in ["Ms", "Mu", "velocity_s_true", "velocity_u_true"]:
        arrays[key] = np.asarray(original.layers[key])
    for prefix, obj in [("default", original), ("gate_open", opened)]:
        for key in ["velo_s", "velo_u"]:
            arrays[f"{prefix}_{key}"] = np.asarray(obj.layers[key])
    np.savez_compressed(output / "inputs/phase_case.npz", **arrays)

    from summarize_controls import summarize
    summaries, effects = summarize(output / "results")
    write(output / "results/condition_summary.csv", summaries)
    write(output / "results/candidate_gate_effects.csv", effects)

    runs = {}
    for name in ["frozen_output_ablation", "all_gene_paired_phase_refits", "velocity_smoothing_ablation", "scvelo_saved_fit_controls_run5", "candidate_gate_pilot", "candidate_gate_branching", "assignment_preview_v3"]:
        raw = json.loads((root / name / "manifest.json").read_text())
        runs[name] = portable(raw)
    external_outputs = []
    for dataset in ["branching", "all"]:
        for condition in ["refit_default", "on_allow_repression"]:
            path = root / "candidate_gate_branching" / f"{dataset}_{condition}.h5ad"
            external_outputs.append(dict(artifact_id=f"multivelo-diagnostic-{dataset}-{condition.replace('_','-')}", locator=portable(str(path)), sha256=sha(path), size_bytes=path.stat().st_size))
    archive_sources = {portable(str(path)): sha(path) for path in copies.values()}
    archive_sources.update({portable(str(root / name)): sha(root / name) for name in ["all_gene_paired_phase_refits/paired_refit_scores.csv", "velocity_smoothing_ablation/paired_refit_scores.csv", "candidate_gate_branching/paired_refit_scores.csv"]})
    provenance = dict(
        scope="Post-hoc diagnostics; primary benchmark results and rankings are unchanged",
        versions=dict(multivelo="0.1.0", scvelo="0.3.3", python="3.9.24", numpy="2.0.2"),
        source_versions={"multivelo_module_sha256": "ac9a998a9b896ed142f0113b19416ef1e0afcf9cc1978487419ca842d9d34020", "multivelo_url": "https://github.com/welch-lab/MultiVelo", "scvelo_url": "https://github.com/theislab/scvelo"},
        reference_neighbors="evaluation/gene/results/reference_neighbors/{dataset}.npz",
        source_runs=runs, source_csvs=archive_sources, external_native_outputs=external_outputs,
        phase_case=dict(dataset="branching", gene="BoostALL_5", n_cells=2000, color="Known progression bins, not continuous physical time", source_conditions=["refit_default", "on_allow_repression"]),
        exclusions=["Hand-reconstructed MultiVelo velocities failed exact output reproduction and are excluded", "scVelo fit_steady_states=False raises an upstream IndexError in this velocity path; no score is reported"],
        outputs={str(p.relative_to(output)): sha(p) for folder in ["results", "inputs"] for p in (output / folder).iterdir() if p.is_file()},
        code_sha256={p.name: sha(p) for p in output.glob("*.py")})
    (output / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    if args.materials_dir is not None:
        publish_materials(args.materials_dir, output)
    print("Published compact results:", output)


if __name__ == "__main__":
    main()
