#!/usr/bin/env python3
"""Replay only the fixed-path DPT step to complete the saved branch audit.

The original main-path values remain the plotting source. Their replay is
checked against the frozen CSV before missing branch values are exported.
No path selection, orientation fitting, velocity fitting or metrics are run.
"""

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-root", type=Path, required=True)
    parser.add_argument("--audit-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=False)
    sys.path.insert(0, str(args.audit_root))
    from stage_audit import SCENARIOS, prepare, set_seed

    frozen_dir = args.audit_root / "ordering_trigger_audit_20260929_threads1"
    frozen_manifest = json.loads((frozen_dir / "manifest.json").read_text())
    frozen_inputs = {row["scenario"]: row for row in frozen_manifest["inputs"]}
    edges_path = args.audit_root / "stage_audit/path_edges_by_stage.csv"
    edges = pd.read_csv(edges_path)
    frozen = pd.read_csv(frozen_dir / "cell_orders.csv", dtype={"cell": str})
    exported, checks, correlations = [], [], []

    for scenario in ("branching", "all"):
        set_seed(frozen_manifest["seed"])
        adata, _, _ = prepare(args.artifact_root, SCENARIOS[scenario])
        for layer in ("Mu", "Ms"):
            digest = hashlib.sha256(np.asarray(adata.layers[layer]).tobytes()).hexdigest()
            if digest != frozen_inputs[scenario][f"{layer}_sha256"]:
                raise ValueError(f"Frozen preprocessing checksum changed: {scenario}/{layer}")
        path_cell_sets = []

        for path in ("main_path", "branch_1"):
            frame = edges[(edges.scenario == scenario) & (edges.stage == "process_path")
                          & (edges.path == path)]
            if frame.empty or list(frame.target.iloc[:-1]) != list(frame.source.iloc[1:]):
                raise ValueError(f"Missing or disconnected frozen path: {scenario}/{path}")
            nodes = list(frame.source) + [frame.iloc[-1].target]
            sub = adata[adata.obs["bin"].astype(str).isin(nodes), :].copy()
            labels = sub.obs["bin"].astype(str).to_numpy()
            # Exact connectivities-mode root-cell rule in TIvelo correct_path.
            sub.uns["iroot"] = int(np.flatnonzero(labels == nodes[0])[0])
            sc.tl.diffmap(sub)
            sc.tl.dpt(sub)
            dpt = sub.obs["dpt_pseudotime"].to_numpy(dtype=float)
            if not np.isfinite(dpt).all():
                raise ValueError(f"Nonfinite path DPT: {scenario}/{path}")
            cells = sub.obs_names.astype(str)
            bins = np.array([int(label.rsplit("_", 1)[1]) for label in labels])
            branches = np.array([label.rsplit("_", 1)[0] for label in labels])
            path_cell_sets.append(set(cells))
            check = {"scenario": scenario, "path": path, "nodes": nodes,
                     "root_cell": str(cells[sub.uns["iroot"]]), "n_cells": len(cells)}

            if path == "main_path":
                old = frozen[(frozen.scenario == scenario) & (frozen.path == "process_main")]
                if set(old.cell) != set(cells):
                    raise ValueError(f"Frozen main-path cell set changed: {scenario}")
                replay = pd.Series(dpt, index=cells)
                maximum_error = float(np.max(np.abs(replay.loc[old.cell].to_numpy() - old.dpt.to_numpy())))
                if maximum_error > 1e-5:
                    raise ValueError(f"Main-path DPT replay differs: {scenario}, max error {maximum_error}")
                check["maximum_main_dpt_replay_error"] = maximum_error
            else:
                exported.extend({"scenario": scenario, "path": path, "cell": cell,
                                 "branch": branch, "bin": int(bin_index), "dpt": float(value)}
                                for cell, branch, bin_index, value in zip(cells, branches, bins, dpt))
                for branch in np.unique(branches):
                    mask = branches == branch
                    correlations.append({"scenario": scenario, "path": path, "branch": branch,
                                         "n_cells": int(mask.sum()),
                                         "spearman_dpt_bin": float(spearmanr(dpt[mask], bins[mask]).statistic)})
            checks.append(check)
            print(json.dumps(check), flush=True)

        if path_cell_sets[0] | path_cell_sets[1] != set(adata.obs_names.astype(str)):
            raise ValueError(f"Saved main and branch paths do not cover all cells: {scenario}")

    pd.DataFrame(exported).to_csv(args.output_dir / "branch_cell_orders.csv", index=False)
    pd.DataFrame(correlations).to_csv(args.output_dir / "branch_dpt_correlations.csv", index=False)
    manifest = {"purpose": "Fixed-path branch DPT export only; no model/metric changes",
                "seed": frozen_manifest["seed"], "source_main_csv": str((frozen_dir / "cell_orders.csv").relative_to(args.audit_root)),
                "path_table_sha256": hashlib.sha256(edges_path.read_bytes()).hexdigest(),
                "inputs": [frozen_inputs[name] for name in ("branching", "all")], "checks": checks,
                "root_rule": "First cell in the initial path node, as correct_path(root_select='connectivities')",
                "independent_path_times": True,
                "generator_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (args.output_dir / "branch_dpt_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(pd.DataFrame(correlations).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
