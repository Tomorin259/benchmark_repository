"""Compare alpha and gamma recovery on the same finite genes in both methods."""

import argparse
from collections import defaultdict
import csv
import hashlib
from io import StringIO
import json
import math
from pathlib import Path
from statistics import fmean
from tempfile import TemporaryDirectory


ROOT = Path(__file__).resolve().parents[2]
SOURCE = Path("evaluation/kinetics_recovery/results/central_results/gene_metrics.csv")
METHODS = ("cellDancer", "DeepVelo")
PARAMETERS = ("alpha", "gamma")
BRANCHES = {"rpe1": ("",), "organoid": ("1", "2")}
SCORE_COLUMNS = ("system", "branch", "gene", "method", "alpha_spearman", "gamma_spearman")


def compare(rows):
    """Select genes with all four finite scores and average within each branch."""
    fits = {(method, parameter) for method in METHODS for parameter in PARAMETERS}
    candidates = defaultdict(dict)
    for row in rows:
        system, branch, gene = row["Dataset"], row["Branch"], row["Gene"]
        fit = (row["Method"], row["Parameter"])
        if system not in BRANCHES or fit not in fits:
            continue
        if branch not in BRANCHES[system] or not gene:
            raise ValueError(f"Invalid gene or branch for {system}: {branch!r}, {gene!r}")
        scores = candidates[system, branch].setdefault(gene, {})
        if fit in scores:
            raise ValueError(f"Duplicate score for {system}/{branch}/{gene}/{fit}")
        scores[fit] = float(row["Spearman_Corr"] or "nan")

    matched_rows = []
    results = {}
    for system, branches in BRANCHES.items():
        components = []
        for branch in branches:
            matched = {
                gene: scores for gene, scores in sorted(candidates[system, branch].items())
                if scores.keys() == fits and all(map(math.isfinite, scores.values()))
            }
            if not matched:
                raise ValueError(f"No jointly finite scores for {system}/{branch or 'all'}")
            matched_rows.extend(
                (system, branch, gene, method, *(scores[method, parameter] for parameter in PARAMETERS))
                for gene, scores in matched.items() for method in METHODS
            )
            means = {
                method: {
                    parameter: fmean(scores[method, parameter] for scores in matched.values())
                    for parameter in PARAMETERS
                }
                for method in METHODS
            }
            components.append({"branch": branch, "common_genes": len(matched), "means": means})

        means = {
            method: {
                parameter: fmean(component["means"][method][parameter] for component in components)
                for parameter in PARAMETERS
            }
            for method in METHODS
        }
        # With two methods, each method is nondominated unless its rival
        # is at least as high on both parameters and strictly higher on one.
        front = []
        for method, rival in (METHODS, METHODS[::-1]):
            differences = [means[rival][parameter] - means[method][parameter]
                           for parameter in PARAMETERS]
            if min(differences) < 0 or max(differences) == 0:
                front.append(method)
        results[system] = {"components": components, "means": means, "front": front}
    return matched_rows, results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="New directory for the comparison.")
    output = parser.parse_args().output
    if output.exists() or output.is_symlink():
        parser.error(f"output directory must be new: {output}")

    source = (ROOT / SOURCE).read_bytes()
    reader = csv.DictReader(StringIO(source.decode("utf-8"), newline=""))
    columns = {"Dataset", "Branch", "Gene", "Method", "Parameter", "Spearman_Corr"}
    missing = columns.difference(reader.fieldnames or ())
    if missing:
        raise ValueError(f"Missing required columns in {SOURCE}: {sorted(missing)}")
    matched_rows, results = compare(reader)
    comparison = {
        "source_file": SOURCE.as_posix(),
        "source_sha256": hashlib.sha256(source).hexdigest(),
        "scope": "Within each system or branch, require finite Spearman scores for alpha and gamma in both methods.",
        "aggregation": "Average genes equally within each component and the two organoid branch means equally.",
        "comparison": "Descriptive Pareto dominance on system-mean alpha and gamma Spearman correlations.",
        "results": results,
    }

    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=f".{output.name}.", dir=output.parent) as temporary:
        staged = Path(temporary) / "comparison"
        staged.mkdir()
        with (staged / "common_gene_scores.csv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, lineterminator="\n")
            writer.writerow(SCORE_COLUMNS)
            writer.writerows(matched_rows)
        (staged / "common_gene_comparison.json").write_text(
            json.dumps(comparison, indent=2, allow_nan=False) + "\n", encoding="utf-8"
        )

        if output.exists() or output.is_symlink():
            parser.error(f"output directory must be new: {output}")
        staged.rename(output)


if __name__ == "__main__":
    main()
