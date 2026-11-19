"""Reproduce supplementary Tables S8E and S10C from committed gene scores."""

import argparse
import csv
from itertools import permutations, product
from pathlib import Path
from statistics import mean

import numpy as np
import pandas as pd

from config import METHOD_ORDER, SIM_CONFIGS

GENE = Path(__file__).resolve().parent
ROOT = GENE.parents[1]
SCENARIOS = ("murk", "branching", "all")
DRAWS = 5_000_000
BATCH = 50_000
SEED = 20261114


def front(scores):
    other, target = scores[:, None, :], scores[None, :, :]
    dominates = (other >= target).all(axis=2) & (other > target).any(axis=2)
    return ~dominates.any(axis=0)


def write_csv(path, rows):
    with path.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def compare(output):
    # Use canonical CBDir and coordinate-aligned CosSim, joined on the same genes.
    reg = pd.read_csv(ROOT / "methods/configurations.tsv", sep="\t")
    labels = {name.casefold(): name for name in reg.manuscript_name}
    source_names = {labels["multivelo (rna_only)" if m == "MultiVelo" else m.casefold()]: m
                    for m in METHOD_ORDER}
    methods = [source_names[name] for name in reg.manuscript_name if name in source_names]
    shown = [name.replace("RNA_only", "RNA only") for name in reg.manuscript_name if name in source_names]
    results = GENE / "results"
    cb = pd.read_csv(results / "gene_cbdir_reference30_all.csv", float_precision="round_trip")
    cb = cb[cb.dataset.isin(SCENARIOS)]
    cos = pd.read_csv(results / "simulation/gene_cbdir_cossim_reference.csv", float_precision="round_trip")
    cos = cos.rename(columns={"Dataset": "dataset", "Method": "method", "Gene": "gene"})
    keys = ["dataset", "method", "gene"]
    scores = cb.merge(cos[keys + ["CosSim"]], on=keys, how="outer", validate="one_to_one")
    if len(scores) != 630 or not np.isfinite(scores[["new_cbdir", "CosSim"]]).all().all() or not scores.status.eq("ok").all():
        raise ValueError("Expected 630 evaluable, paired simulation scores")
    scenarios, subgroup = [], {}
    for dataset in SCENARIOS:
        part = scores[scores.dataset == dataset]
        expected = set(SIM_CONFIGS[dataset]["genes"])
        if any(set(part.loc[part.method == method, "gene"]) != expected for method in methods):
            raise ValueError(f"Incomplete common gene set in {dataset}")
        scenarios.append(part.groupby("method")[["new_cbdir", "CosSim"]].agg(
            lambda values: mean(map(float, values))).loc[methods].to_numpy())
        if dataset == "all":
            gene_sets = {"murk": part.gene.str.startswith(("OneSidedBoost_", "OneSidedBoost2_")),
                         "branching": part.gene.str.startswith(("MultiForward_", "BoostALL_"))}
            for kind, mask in gene_sets.items():
                if part.loc[mask, "gene"].nunique() != 20:
                    raise ValueError(f"Expected 20 hybrid {kind} genes")
                subgroup[kind] = part[mask].groupby("method").new_cbdir.mean().loc[methods].to_numpy()
    cube = np.stack(scenarios, axis=1)  # method, scenario, metric
    membership = front(cube.reshape(len(methods), 6))
    draws = list(product(range(3), repeat=3))
    hits = sum(front(cube[:, draw, :].reshape(len(methods), 6)) for draw in draws)
    s8e = []
    for i, method in enumerate(shown):
        row = {"method": method}
        for j, prefix in enumerate(("normal_murk", "normal_branching", "hybrid")):
            row[prefix + "_gene_cbdir"], row[prefix + "_cossim"] = cube[i, j]
        row.update(on_full_front=membership[i], front_hits=hits[i], ordered_draws=len(draws))
        s8e.append(row)

    real = pd.read_csv(results / "common_genes/dataset_summary.csv")
    real = real.pivot(index="dataset", columns="method", values="mean")[methods]
    selected = pd.read_csv(GENE / "selection/selected_genes.tsv", sep="\t")
    kinds = selected[["dataset", "setting"]].drop_duplicates().set_index("dataset").setting
    if not kinds.index.is_unique or set(kinds.index) != set(real.index) or real.shape != (7, 9):
        raise ValueError("Expected the registered seven real datasets and nine methods")
    murk = real.loc[kinds[kinds == "murk"].index].to_numpy()
    branching = real.loc[kinds[kinds == "branching"].index].to_numpy()
    groups = [
        ("source_real_cbdir", "real", "CBDir", real.to_numpy()),
        ("source_sim_cbdir", "simulation", "CBDir", cube[:, :, 0].T),
        ("source_sim_cossim", "simulation", "CosSim", cube[:, :, 1].T),
        ("category_murk_cbdir", "MURK", "CBDir", np.vstack((murk, cube[:, 0, 0], subgroup["murk"]))),
        ("category_branching_cbdir", "branching", "CBDir", np.vstack((branching, cube[:, 1, 0], subgroup["branching"]))),
    ]
    tests, thresholds = [], []
    for name, group, metric, matrix in groups:
        n, k = matrix.shape
        if not np.isfinite(matrix).all() or any(len(set(row)) != k for row in matrix):
            raise ValueError(f"The frozen untied-rank permutation model does not apply to {name}")
        ranks = (-matrix).argsort(axis=1).argsort(axis=1) + 1
        threshold = int(np.square((2 * ranks - (k + 1)).sum(axis=0)).sum())
        q = 3 * threshold / (n * k * (k + 1))
        thresholds.append(threshold)
        tests.append(dict(analysis_id=name, group=group, metric=metric, n_blocks=n, n_methods=k,
                          friedman_q=q, kendall_w=q / (n * (k - 1)), permutations=DRAWS,
                          extreme_permutations=0, seed=SEED))

    # Relabeling all columns together leaves Q unchanged, so fix the first block.
    base = np.arange(-8, 9, 2, dtype=np.int16)
    pool = np.asarray(list(permutations(base)), dtype=np.int16)
    for n in sorted({row["n_blocks"] for row in tests}):
        family = [i for i, row in enumerate(tests) if row["n_blocks"] == n]
        rng = np.random.default_rng(SEED)
        for start in range(0, DRAWS, BATCH):
            size = min(BATCH, DRAWS - start)
            totals = np.broadcast_to(base, (size, len(base))).copy()
            for _ in range(n - 1):
                totals += pool[rng.integers(0, len(pool), size=size)]
            statistic = np.square(totals.astype(np.int64)).sum(axis=1)
            for i in family:
                tests[i]["extreme_permutations"] += int((statistic >= thresholds[i]).sum())
    for row in tests:
        row["p_mc"] = (row["extreme_permutations"] + 1) / (DRAWS + 1)
    previous = 0.0
    for index, row in enumerate(sorted(tests, key=lambda r: r["p_mc"])):
        previous = max(previous, min(1.0, (len(tests) - index) * row["p_mc"]))
        row["p_holm_global5"] = previous
    output.mkdir(parents=True, exist_ok=False)
    write_csv(output / "s8e_six_objective_front.csv", s8e)
    write_csv(output / "s10c.csv", [row for row in tests if row["group"] == "simulation"])
    write_csv(output / "holm5_calibration.csv", tests)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True, help="A new directory for the three result tables")
    compare(parser.parse_args().output_dir)
