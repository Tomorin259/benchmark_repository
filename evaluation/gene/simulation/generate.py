"""Generate additional datasets for the three study simulation scenarios.

Scenario parameters and randomness settings are recorded in config.json.
"""
import argparse
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
CONFIG_PATH = HERE / "config.json"


def load_config():
    return json.loads(CONFIG_PATH.read_text())


def _sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _simulate(cs, kinetic_type, kwargs):
    import numpy as np
    original_jitter = cs._jitter
    clean = []

    def retain_state(u, s, scale):
        clean.append((np.asarray(u).copy(), np.asarray(s).copy()))
        return original_jitter(u, s, scale)

    cs._jitter = retain_state
    try:
        frame = cs.simulate(kinetic_type=kinetic_type, **kwargs)
    finally:
        cs._jitter = original_jitter
    if not clean or sum(len(u) for u, _ in clean) != len(frame):
        raise ValueError("cellDancer noise call does not match generated cells")
    frame["clean_unsplice"] = np.concatenate([u for u, _ in clean])
    frame["clean_splice"] = np.concatenate([s for _, s in clean])
    return frame


def _base_kwargs(scenario):
    return {
        "path1_cell_number": scenario["path_cells"][0],
        "path2_cell_number": scenario["path_cells"][1],
        "path1_pct": scenario["path_percent"][0],
        "path2_pct": scenario["path_percent"][1],
        "noise_level": scenario["noise_level"],
    }


def _make_gene(cs, scenario, parameters):
    common = _base_kwargs(scenario)
    return _simulate(cs, parameters["kinetic_type"], dict(common, **{k: v for k, v in parameters.items() if k != "kinetic_type"}))


def _one_sided(cs, scenario, parameters, ranges):
    import numpy as np
    common = _base_kwargs(scenario)
    p1, p2 = scenario["path_cells"]
    alpha1 = float(np.random.uniform(*ranges[0]))
    alpha2 = float(np.random.uniform(*ranges[1]))
    quiet = _simulate(cs, parameters["quiet"]["kinetic_type"], dict(common, **{k: v for k, v in parameters["quiet"].items() if k != "kinetic_type"},
        path1_cell_number=p1 // 2, path2_cell_number=p1 - p1 // 2))
    boost = _simulate(cs, parameters["boost"]["kinetic_type"], dict(common, alpha1=alpha1, alpha2=alpha2,
        **{k: v for k, v in parameters["boost"].items() if k not in {"kinetic_type", "path1_pct", "path2_pct"}},
        path1_cell_number=p2 // 2, path2_cell_number=p2 - p2 // 2, path1_pct=99, path2_pct=90))
    import pandas as pd
    frame = pd.concat([quiet, boost], ignore_index=True)
    frame["cellID"] = np.arange(p1 + p2)
    return frame


def _bin_labels(cell_ids, path_cells, retain):
    if retain == "branch_b":
        return ["BranchB_%d" % (int(cell_id) // 100) for cell_id in cell_ids], ["BranchB_%d" % i for i in range(10)]
    labels = []
    for cell_id in cell_ids:
        cell_id = int(cell_id)
        labels.append("BranchA_%d" % (cell_id // 100) if cell_id < path_cells[0] else "BranchB_%d" % ((cell_id - path_cells[0]) // 100))
    return labels, ["BranchA_%d" % i for i in range(10)] + ["BranchB_%d" % i for i in range(10)]


def _truth_layers(s, u, alpha, beta, gamma):
    return beta * u - gamma * s, alpha - beta * u


def _runtime_versions():
    import pkg_resources
    return {
        name: pkg_resources.get_distribution(name).version
        for name in ("numpy", "scipy", "pandas", "anndata", "celldancer")
    }


def generate(scenario_name, output_root):
    import numpy as np
    config = load_config()
    if scenario_name not in config["scenarios"]:
        raise ValueError("scenario must be one of: " + ", ".join(config["scenarios"]))
    scenario = config["scenarios"][scenario_name]
    parameters = config["gene_parameters"]
    output_root = Path(output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    import anndata
    import celldancer.simulation as cs
    import pandas as pd

    # CellDancer 1.1.7 and the parameter sampler use NumPy's global random stream.
    frames = []
    counts = scenario["genes"]
    stem = None
    for kind in ("housekeeping", "boost_all", "one_sided_boost", "one_sided_boost2", "stem_marker", "multi_forward"):
        for i in range(counts[kind]):
            if kind == "one_sided_boost":
                frame = _one_sided(cs, scenario, {"quiet": parameters["one_sided_quiet"], "boost": parameters["one_sided_boost"]}, scenario["one_sided_alpha_ranges"])
            elif kind == "one_sided_boost2":
                frame = _one_sided(cs, scenario, {"quiet": parameters["one_sided_quiet"], "boost": parameters["one_sided_boost"]}, scenario["one_sided_alpha_ranges2"])
                frame = pd.concat([frame.iloc[scenario["path_cells"][1]:], frame.iloc[:scenario["path_cells"][1]]], ignore_index=True)
                frame["cellID"] = np.arange(sum(scenario["path_cells"]))
            elif kind == "stem_marker" and stem is not None:
                frame = stem.copy()
            else:
                frame = _make_gene(cs, scenario, parameters[kind])
                if kind == "stem_marker":
                    stem = frame.copy()
            frame["gene_name"] = "%s_%d" % (kind.replace("_", "").title(), i)
            if kind == "boost_all":
                frame["gene_name"] = "BoostALL_%d" % i
            elif kind == "one_sided_boost":
                frame["gene_name"] = "OneSidedBoost_%d" % i
            elif kind == "one_sided_boost2":
                frame["gene_name"] = "OneSidedBoost2_%d" % i
            elif kind == "stem_marker":
                frame["gene_name"] = "StemMarker_%d" % i
            elif kind == "multi_forward":
                frame["gene_name"] = "MultiForward_%d" % i
            frames.append(frame)

    full = pd.concat(frames, ignore_index=True)
    if scenario["retain"] == "branch_b":
        full = full[full["cellID"] >= scenario["path_cells"][0]].copy()
        full["cellID"] -= scenario["path_cells"][0]
    genes = sorted(full["gene_name"].unique())
    s = full.pivot(index="cellID", columns="gene_name", values="splice").reindex(columns=genes)
    u = full.pivot(index="cellID", columns="gene_name", values="unsplice").reindex(columns=genes)
    clean_s = full.pivot(index="cellID", columns="gene_name", values="clean_splice").reindex(columns=genes)
    clean_u = full.pivot(index="cellID", columns="gene_name", values="clean_unsplice").reindex(columns=genes)
    alpha = full.pivot(index="cellID", columns="gene_name", values="alpha").reindex(columns=genes)
    beta = full.pivot(index="cellID", columns="gene_name", values="beta").reindex(columns=genes)
    gamma = full.pivot(index="cellID", columns="gene_name", values="gamma").reindex(columns=genes)
    obs = pd.DataFrame(index=[str(x) for x in s.index])
    obs["true_branch"] = "Branch_B" if scenario["retain"] == "branch_b" else ["Branch_A" if int(x) < scenario["path_cells"][0] else "Branch_B" for x in s.index]
    labels, categories = _bin_labels(s.index, scenario["path_cells"], scenario["retain"])
    obs["bin"] = pd.Categorical(labels, categories=categories, ordered=True)
    obs["clusters"] = obs["true_branch"].copy()
    adata = anndata.AnnData(X=s.to_numpy(), obs=obs, var=pd.DataFrame(index=genes))
    adata.layers["spliced"] = s.to_numpy()
    adata.layers["unspliced"] = u.to_numpy()
    adata.layers["alpha_true"] = alpha.to_numpy()
    adata.layers["beta_true"] = beta.to_numpy()
    adata.layers["gamma_true"] = gamma.to_numpy()
    adata.layers["clean_spliced"] = clean_s.to_numpy()
    adata.layers["clean_unspliced"] = clean_u.to_numpy()
    adata.layers["velocity_s_true"], adata.layers["velocity_u_true"] = _truth_layers(clean_s.to_numpy(), clean_u.to_numpy(), alpha.to_numpy(), beta.to_numpy(), gamma.to_numpy())
    if list(adata.shape) != scenario["expected_shape"]:
        raise ValueError("generated shape %s does not match expected %s" % (adata.shape, scenario["expected_shape"]))
    output = output_root / scenario["output_name"]
    adata.write_h5ad(output)
    metadata = {
        "schema_version": 1, "scenario": scenario_name,
        "rng_initialization": "numpy.random global RandomState; no explicit seed",
        "config_sha256": _sha256(CONFIG_PATH),
        "generator_sha256": _sha256(__file__),
        "output": output.name, "shape": list(adata.shape),
        "output_sha256": _sha256(output),
        "runtime_versions": _runtime_versions(),
        "canonical_artifact_id": scenario["artifact_id"],
        "velocity_reference": "ODE velocity at the retained pre-noise state",
        "realization": "independent"
    }
    (output.with_suffix(".invocation.json")).write_text(json.dumps(metadata, indent=2) + "\n")
    return output


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--scenario", required=True, choices=sorted(load_config()["scenarios"]))
    parser.add_argument("--output-root", default="outputs/simulation")
    args = parser.parse_args()
    generate(args.scenario, args.output_root)


if __name__ == "__main__":
    main()
