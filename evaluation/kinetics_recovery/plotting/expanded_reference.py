#!/usr/bin/env python3
"""Render expanded-reference materials through the shared Figure 6B/6F functions."""
from __future__ import annotations

import argparse
from pathlib import Path
import tempfile

from ..expanded_reference import CANONICAL, candidate_genes, load_scores
from .figure6 import (DISPLAY_METHOD_NAMES, canonical_method_color, ensure_output,
                      plot_branch, plot_trajectory)


def render(output: Path, *, trajectories: bool = True) -> None:
    scope, scores = load_scores()
    with tempfile.TemporaryDirectory(prefix="expanded-kinetics-") as directory:
        score_path = Path(directory) / "scores.csv"
        scores.to_csv(score_path, index=False)
        for dataset, branch, label, filename in (
            ("rpe1", None, "RPE1", "rpe1_boxplot.png"),
            ("organoid", 1, "Organoid - Branch 1", "organoid_branch1_boxplot.png"),
            ("organoid", 2, "Organoid - Branch 2", "organoid_branch2_boxplot.png"),
        ):
            frame = scores[scores.Dataset.eq(dataset)]
            frame = frame[frame.Branch.isna()] if branch is None else frame[frame.Branch.eq(branch)]
            paths = {}
            for method in ("Dynamo", "DeepKINET", "TSvelo", "cellDancer", "DeepVelo"):
                group = frame[frame.Method.eq(DISPLAY_METHOD_NAMES[method])]
                if group.empty:
                    continue
                paths[method] = {"color": canonical_method_color(method), **{
                    parameter.capitalize(): score_path for parameter in group.Parameter.unique()
                }}
            plotted = plot_branch(
                branch, paths, candidate_genes(scope, dataset), output / filename,
                dataset=dataset, system_label=label, show_n=True, point_size=2.5,
                point_alpha=0.25, point_color="#737373", point_linewidth=0,
                dpi=600, bbox_inches=None,
            )
            # Exact plotted rows are reconstructible; retain them only in new rendering runs.
            plotted.to_csv(output / f"{Path(filename).stem}_plot_data.csv", index=False)
        if trajectories:
            from matplotlib.backends.backend_pdf import PdfPages
            genes = candidate_genes(scope, "rpe1")
            with PdfPages(output / "rpe1_250gene_trajectories.pdf") as writer:
                for parameter in ("alpha", "gamma"):
                    frame = scores[scores.Dataset.eq("rpe1") & scores.Parameter.eq(parameter)]
                    methods = {
                        method: CANONICAL / f"trajectories/rpe1/{method}_{parameter}.csv"
                        for method, display in DISPLAY_METHOD_NAMES.items() if display in set(frame.Method)
                    }
                    for first in range(0, len(genes), 50):
                        plot_trajectory(
                            parameter, CANONICAL / f"matrices/rpe1/reference_{parameter}.csv",
                            methods, genes[first:first+50], writer, complete_legend=True,
                        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--skip-trajectories", action="store_true")
    args = parser.parse_args()
    ensure_output(args.output)
    render(args.output, trajectories=not args.skip_trajectories)


if __name__ == "__main__":
    main()
