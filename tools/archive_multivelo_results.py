"""Import a complete, verified MultiVelo diagnostic bundle or its materials."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import shutil

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def publish_materials(source, output):
    record = json.loads((source / "manifest.json").read_text())
    if record["script_sha256"] != sha(output / "plot_materials.py"):
        raise ValueError("Material renderer differs from the current source")
    if not record["input_hashes_unchanged"]:
        raise ValueError("Material inputs changed during rendering")
    for name, expected in record["inputs"].items():
        if sha(output / name) != expected:
            raise ValueError(f"Material input differs: {name}")
    published = {}
    for name, expected in record["outputs"].items():
        if Path(name).name != name or sha(source / name) != expected:
            raise ValueError(f"Invalid material identity: {name}")
        folder = "tables" if name == "tableS9.tex" else "figures"
        target = output / folder / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source / name, target)
        published[str(target.relative_to(output))] = sha(target)
    record["outputs"] = published
    record["scope"] = "Source materials for Supplementary Figure S15 A/B and Table S9"
    (output / "materials_manifest.json").write_text(json.dumps(record, indent=2) + "\n")


def import_results(source, output):
    record = json.loads((source / "provenance.json").read_text())
    reference = ROOT / "evaluation/gene/simulation/truth/manifest.json"
    if record.get("velocity_reference", {}).get("sha256") != sha(reference):
        raise ValueError("Diagnostic bundle does not use the current ODE reference")
    required = {"inputs/phase_case.npz", "results/scvelo_controls.csv", "results/multivelo_refits.csv",
                "results/condition_summary.csv", "results/candidate_gate_branch_scores.csv",
                "results/candidate_gate_pilot_cell_states.csv", "results/constant_direction_baselines.csv"}
    if not required.issubset(record.get("outputs", {})):
        raise ValueError("Incomplete diagnostic bundle output manifest")
    for name, expected_rows, keys in [
        ("scvelo_controls.csv", 490, ["dataset", "gene", "variant"]),
        ("multivelo_refits.csv", 250, ["dataset", "gene", "condition"]),
        ("condition_summary.csv", 22, ["dataset", "gene_class", "method", "condition"]),
        ("candidate_gate_branch_scores.csv", 160, ["dataset", "gene", "condition", "branch"]),
        ("candidate_gate_pilot_cell_states.csv", 4000, ["condition", "cell_id"]),
        ("constant_direction_baselines.csv", 210, ["dataset", "gene", "baseline"]),
    ]:
        with (source / "results" / name).open(newline="") as handle:
            rows = list(csv.DictReader(handle))
        if len(rows) != expected_rows or len({tuple(row[key] for key in keys) for row in rows}) != expected_rows:
            raise ValueError(f"Incomplete or duplicated diagnostic comparison: {name}")
    for name, expected in record["outputs"].items():
        if Path(name).parts[0] not in {"results", "inputs"} or ".." in Path(name).parts:
            raise ValueError(f"Invalid result path: {name}")
        if sha(source / name) != expected:
            raise ValueError(f"Diagnostic output differs: {name}")
    for name in record["outputs"]:
        destination = output / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source / name, destination)
    (output / "provenance.json").write_text(json.dumps(record, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--diagnostics-root", type=Path, help="Complete bundle with results/, inputs/ and provenance.json")
    parser.add_argument("--output-root", type=Path, default=ROOT / "supplement/multivelo")
    parser.add_argument("--materials-dir", type=Path)
    parser.add_argument("--materials-only", action="store_true")
    args = parser.parse_args()
    if args.materials_only and (args.materials_dir is None or args.diagnostics_root is not None):
        parser.error("--materials-only requires --materials-dir and excludes --diagnostics-root")
    if not args.materials_only and args.diagnostics_root is None:
        parser.error("Supply a complete --diagnostics-root bundle")
    output = args.output_root.resolve()
    output.mkdir(parents=True, exist_ok=True)
    if not args.materials_only:
        import_results(args.diagnostics_root.resolve(), output)
    if args.materials_dir:
        publish_materials(args.materials_dir.resolve(), output)
    print("Imported verified diagnostics:", output)


if __name__ == "__main__":
    main()
