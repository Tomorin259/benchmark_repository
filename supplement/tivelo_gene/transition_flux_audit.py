#!/usr/bin/env python3
"""Summarize scVelo transition flux along the simulator's annotated bins.

The saved TIvelo H5AD contains the upstream stochastic scVelo velocity graph;
this script only reads it and constructs scVelo's forward transition matrix.
It does not rerun or modify the benchmark.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scvelo.tools.transition_matrix import transition_matrix


LABEL = re.compile(r"^Branch([AB])_(\d+)$")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent / "transition_flux_audit")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows, pair_rows = [], []
    for scenario in ("murk", "branching", "all"):
        path = args.artifact_root / f"indicators/gene/simulation/results_for_analyze/tivelo/tivelo_{scenario}.h5ad"
        adata = ad.read_h5ad(path)
        T = transition_matrix(adata, vkey="velocity", backward=False, self_transitions=False).tocsr()
        labels = adata.obs["bin"].astype(str).to_numpy()
        categories = list(adata.obs["bin"].cat.categories)
        for source in categories:
            source_rows = np.flatnonzero(labels == source)
            lmatch = LABEL.fullmatch(source)
            if lmatch is None:
                continue
            branch, time = lmatch.group(1), int(lmatch.group(2))
            mass = T[source_rows, :].sum(axis=0).A.ravel() / len(source_rows)
            total_forward = total_reverse = total_same = total_cross = 0.0
            for target in categories:
                rmatch = LABEL.fullmatch(target)
                target_mass = float(mass[labels == target].sum())
                pair_rows.append({"scenario": scenario, "source": source, "target": target,
                                  "mean_transition_probability": target_mass})
                if rmatch.group(1) != branch:
                    total_cross += target_mass
                elif int(rmatch.group(2)) > time:
                    total_forward += target_mass
                elif int(rmatch.group(2)) < time:
                    total_reverse += target_mass
                else:
                    total_same += target_mass
            rows.append({"scenario": scenario, "source": source, "n_cells": len(source_rows),
                         "forward_probability": total_forward,
                         "reverse_probability": total_reverse,
                         "same_bin_probability": total_same,
                         "cross_branch_probability": total_cross,
                         "forward_minus_reverse": total_forward - total_reverse,
                         "mean_root_score": float(np.mean(adata.obs["root_cells"].to_numpy()[source_rows])),
                         "mean_end_score": float(np.mean(adata.obs["end_points"].to_numpy()[source_rows]))})
        print(f"completed {scenario}", flush=True)
    pd.DataFrame(rows).to_csv(args.output_dir / "transition_flux_by_source_bin.csv", index=False)
    pd.DataFrame(pair_rows).to_csv(args.output_dir / "transition_probability_by_bin_pair.csv", index=False)
    (args.output_dir / "manifest.json").write_text(json.dumps({
        "purpose": "scVelo transition flux upstream of TIvelo root selection",
        "inputs": [f"artifact://velocity/indicators/gene/simulation/results_for_analyze/tivelo/tivelo_{scenario}.h5ad"
                   for scenario in ("murk", "branching", "all")],
        "transition_matrix": "scvelo.tools.transition_matrix(vkey='velocity', backward=False, self_transitions=False)",
    }, indent=2) + "\n")
    print(pd.DataFrame(rows).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
