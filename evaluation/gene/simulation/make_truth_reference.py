"""Materialize noise-free ODE states and velocities for registered Gene cases."""

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def simulate(cs, kinetic_type, rates, n1, n2, pct1, pct2):
    alpha, beta, gamma = rates
    return cs.simulate(
        kinetic_type=kinetic_type,
        alpha1=float(alpha[0]), alpha2=float(alpha[n1]),
        beta1=float(beta[0]), beta2=float(beta[n1]),
        gamma1=float(gamma[0]), gamma2=float(gamma[n1]),
        path1_cell_number=n1, path2_cell_number=n2,
        path1_pct=pct1, path2_pct=pct2, noise_level=0,
    )


def make_reference(repository_root, artifact_root, output_root, all_genes=False):
    import anndata as ad
    import celldancer.simulation as cs
    import pandas as pd
    import pkg_resources

    sys.path.insert(0, str(repository_root))
    from evaluation.gene.config import SIM_CONFIGS
    from evaluation.gene.simulation.generate import _simulate

    config = json.loads((repository_root / "evaluation/gene/simulation/config.json").read_text())
    with (repository_root / "artifacts.tsv").open(newline="") as handle:
        registered = {row["artifact_id"]: row["sha256"] for row in csv.DictReader(handle, delimiter="\t")}
    if output_root.exists() and any(output_root.iterdir()):
        raise ValueError("Truth output directory must be empty")
    output_root.mkdir(parents=True, exist_ok=True)
    records = []
    for scenario in ("murk", "branching", "all"):
        source = artifact_root / SIM_CONFIGS[scenario]["reference_source"]
        config_name = "hybrid" if scenario == "all" else scenario
        artifact_id = config["scenarios"][config_name]["artifact_id"]
        source_sha = sha256(source)
        if source_sha != registered[artifact_id]:
            raise ValueError("Registered simulation input hash differs for " + scenario)
        reference = ad.read_h5ad(source, backed="r")
        try:
            cells = np.asarray(reference.obs_names).astype(str)
            genes = np.asarray(reference.var_names if all_genes else SIM_CONFIGS[scenario]["genes"]).astype(str)
            if len(set(cells)) != len(cells) or len(set(genes)) != len(genes):
                raise ValueError("Duplicate cell or gene identifier in " + scenario)
            if cells.tolist() != [str(i) for i in range(len(cells))]:
                raise ValueError("Unexpected cell order in " + scenario)
            positions = reference.var_names.get_indexer(genes)
            if np.any(positions < 0):
                raise ValueError("Missing selected gene in " + scenario)
            clean_u = np.empty((len(cells), len(genes)), dtype=np.float64)
            clean_s = np.empty_like(clean_u)
            velocity_u = np.empty_like(clean_u)
            velocity_s = np.empty_like(clean_u)
            for column, gene in enumerate(genes):
                position = positions[column]
                rates = tuple(
                    np.asarray(reference.layers[key][:, position], dtype=np.float64).ravel()
                    for key in ("alpha_true", "beta_true", "gamma_true")
                )
                alpha, beta, gamma = rates
                family = gene.rsplit("_", 1)[0]
                if family in ("Housekeeping", "StemMarker"):
                    settings = config["scenarios"][config_name]
                    n1, n2 = settings["path_cells"]
                    parameters = config["gene_parameters"]["housekeeping" if family == "Housekeeping" else "stem_marker"]
                    simulated = _simulate(cs, parameters["kinetic_type"], dict(
                        **{key: value for key, value in parameters.items() if key != "kinetic_type"},
                        path1_cell_number=n1, path2_cell_number=n2,
                        path1_pct=settings["path_percent"][0], path2_pct=settings["path_percent"][1], noise_level=0))
                    simulated["splice"] = simulated["clean_splice"]
                    simulated["unsplice"] = simulated["clean_unsplice"]
                    if scenario == "murk":
                        simulated = simulated.iloc[n1:].reset_index(drop=True)
                elif family in ("MultiForward", "BoostALL"):
                    simulated = simulate(cs, "multi_forward", rates, 1000, 1000, 95, 95)
                elif family in ("OneSidedBoost", "OneSidedBoost2"):
                    if scenario == "murk":
                        simulated = simulate(cs, "tran_boost", rates, 500, 500, 99, 90)
                    else:
                        start = 1000 if family == "OneSidedBoost" else 0
                        boost_rates = tuple(rate[start:start + 1000] for rate in rates)
                        boost = simulate(cs, "tran_boost", boost_rates, 500, 500, 99, 90)
                        quiet = cs.simulate(
                            kinetic_type="mono", alpha1=0.5, alpha2=0.5,
                            beta1=1, beta2=1, gamma1=0.5, gamma2=0.5,
                            path1_cell_number=500, path2_cell_number=500,
                            path1_pct=99, path2_pct=99, noise_level=0,
                        )
                        parts = [quiet, boost] if family == "OneSidedBoost" else [boost, quiet]
                        simulated = pd.concat(parts, ignore_index=True)
                else:
                    raise ValueError("Unknown selected gene family: " + family)
                if len(simulated) != len(cells):
                    raise ValueError("Cell count differs for " + scenario + "/" + gene)
                for field, expected in zip(("alpha", "beta", "gamma"), rates):
                    if not np.allclose(simulated[field].to_numpy(), expected, rtol=0, atol=1e-12):
                        raise ValueError("Kinetic parameter alignment differs for " + scenario + "/" + gene)
                clean_u[:, column] = simulated["unsplice"].to_numpy()
                clean_s[:, column] = simulated["splice"].to_numpy()
                velocity_u[:, column] = alpha - beta * clean_u[:, column]
                velocity_s[:, column] = beta * clean_u[:, column] - gamma * clean_s[:, column]
            output = output_root / (scenario + "_truth.npz")
            np.savez_compressed(
                output, cell_ids=cells, gene_names=genes,
                clean_u=clean_u, clean_s=clean_s,
                velocity_u_true=velocity_u, velocity_s_true=velocity_s,
            )
            records.append({
                "scenario": scenario,
                "source": SIM_CONFIGS[scenario]["reference_source"],
                "source_sha256": source_sha,
                "output": output.name,
                "output_sha256": sha256(output),
                "shape": [len(cells), len(genes)],
                "cell_ids_unique": True,
                "kinetic_parameters_equal": True,
            })
            print(scenario, len(cells), len(genes), records[-1]["output_sha256"], flush=True)
        finally:
            reference.file.close()
    manifest = {
        "schema_version": 1,
        "truth": "noise-free cellDancer ODE state",
        "cellDancer_version": "1.1.7",
        "generator_sha256": sha256(Path(__file__)),
        "config_sha256": sha256(repository_root / "evaluation/gene/simulation/config.json"),
        "runtime_versions": {name: pkg_resources.get_distribution(name).version
                             for name in ("celldancer", "numpy", "scipy", "pandas", "anndata")},
        "coordinate_order": ["spliced", "unspliced"],
        "scenario_records": records,
    }
    (output_root / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository-root", type=Path, default=ROOT)
    parser.add_argument("--artifact-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--all-genes", action="store_true", help="Include regular background genes for full-expression diagnostics")
    args = parser.parse_args()
    make_reference(args.repository_root.resolve(), args.artifact_root.resolve(), args.output_root.resolve(), args.all_genes)


if __name__ == "__main__":
    main()
