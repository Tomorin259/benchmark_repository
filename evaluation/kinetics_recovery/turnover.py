#!/usr/bin/env python3
"""Reproduce Supplementary Table S11 from the four Figure 6D CSVs."""

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np


PLOT_DATA = Path(__file__).resolve().parent / "figures/fig6/plot_data"
METHODS = ("cellDancer", "DeepVelo")
CLASSES = ("Cooperative", "Neutral", "Destabilizing")
N_GENES = 250
N_BOOTSTRAP = 10_000
N_PERMUTATIONS = 50_000
SEED = 20261111
SUFFIXES = ("kinetic_class_membership", "kinetic_classes")


def read_csv(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def class_labels(membership_file, counts_file, method):
    reference, prediction = {}, {}
    for row in read_csv(membership_file):
        index = CLASSES.index(row["Class"])
        for column, mapping in (("scEU", reference), (method, prediction)):
            value = row[column]
            gene = row["Gene"]
            if value not in ("True", "False") or not gene:
                raise ValueError(f"{membership_file}: invalid membership for {column}, {gene!r}")
            if value == "True":
                if gene in mapping:
                    raise ValueError(f"{membership_file}: duplicate {column} label for {gene}")
                mapping[gene] = index
    if set(reference) != set(prediction) or len(reference) != N_GENES:
        raise ValueError(f"{membership_file}: expected {N_GENES} matched genes")

    counts = read_csv(counts_file)
    if sorted(row["Class"] for row in counts) != sorted(CLASSES):
        raise ValueError(f"{counts_file}: expected exactly one row per class")
    for row in counts:
        index = CLASSES.index(row["Class"])
        true_set = {g for g, c in reference.items() if c == index}
        predicted_set = {g for g, c in prediction.items() if c == index}
        actual = (
            len(true_set - predicted_set),
            len(predicted_set - true_set),
            len(true_set & predicted_set),
        )
        expected = tuple(int(row[column]) for column in ("scEU_only", f"{method}_only", "intersection"))
        if actual != expected:
            raise ValueError(f"{counts_file}: {row['Class']} counts disagree with {membership_file}")
    return reference, prediction


def class_f1(reference, prediction):
    values = []
    for index in range(len(CLASSES)):
        true_mask = reference == index
        predicted_mask = prediction == index
        positive = np.count_nonzero(true_mask & predicted_mask, axis=-1)
        total = np.count_nonzero(true_mask, axis=-1) + np.count_nonzero(predicted_mask, axis=-1)
        values.append(np.divide(2 * positive, total, out=np.zeros_like(positive, dtype=float), where=total > 0))
    return np.stack(values, axis=-1)


def include_macro(values):
    return np.concatenate([values, values.mean(axis=-1, keepdims=True)], axis=-1)


def turnover(inputs):
    reference, celldancer = class_labels(*inputs[METHODS[0]], METHODS[0])
    other_reference, deepvelo = class_labels(*inputs[METHODS[1]], METHODS[1])
    if reference != other_reference:
        raise ValueError("Fig. 6D scEU-seq reference memberships differ between methods")
    genes = sorted(reference)
    labels = np.array([reference[gene] for gene in genes])
    predictions = [np.array([mapping[gene] for gene in genes]) for mapping in (celldancer, deepvelo)]
    observed = [include_macro(class_f1(labels, p)) for p in predictions]
    rng_boot, rng_perm = [np.random.default_rng(seed) for seed in np.random.SeedSequence(SEED).spawn(2)]
    strata = [np.flatnonzero(labels == i) for i in range(len(CLASSES))]
    draws = np.concatenate([
        rng_boot.choice(group, size=(N_BOOTSTRAP, len(group)), replace=True)
        for group in strata
    ], axis=1)
    intervals = [np.quantile(include_macro(class_f1(labels[draws], p[draws])), [0.025, 0.975], axis=0)
                 for p in predictions]
    diff = observed[1] - observed[0]
    as_extreme = np.zeros(len(CLASSES), dtype=int)
    for start in range(0, N_PERMUTATIONS, 2000):
        batch = min(2000, N_PERMUTATIONS - start)
        swap = rng_perm.integers(0, 2, size=(batch, len(genes)), dtype=np.int8).astype(bool)
        first = np.where(swap, predictions[1], predictions[0])
        second = np.where(swap, predictions[0], predictions[1])
        permuted = class_f1(labels, second) - class_f1(labels, first)
        as_extreme += (np.abs(permuted) >= np.abs(diff[:-1]) - 1e-12).sum(axis=0)
    raw = (as_extreme + 1) / (N_PERMUTATIONS + 1)
    order = np.argsort(raw)
    adjusted = np.empty_like(raw)
    adjusted[order] = np.minimum(1, np.maximum.accumulate(raw[order] * np.arange(len(CLASSES), 0, -1)))
    rows = []
    for i, name in enumerate((*CLASSES, "Macro")):
        row = {
            "class": name,
            "reference_genes": len(strata[i]) if i < len(CLASSES) else len(genes),
        }
        for j, method in enumerate(METHODS):
            row[f"{method}_f1"] = float(observed[j][i])
            row[f"{method}_ci95_low"] = float(intervals[j][0, i])
            row[f"{method}_ci95_high"] = float(intervals[j][1, i])
        row.update({
            "difference_deepvelo_minus_celldancer": float(diff[i]),
            "permutation_p": float(raw[i]) if i < len(CLASSES) else None,
            "holm_p": float(adjusted[i]) if i < len(CLASSES) else None,
        })
        rows.append(row)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=PLOT_DATA,
                        help="Directory containing the four Figure 6D membership and count CSVs")
    parser.add_argument("--output", type=Path, required=True, help="New output directory")
    args = parser.parse_args()

    names = {
        method: [f"rpe1_{method.lower()}_{suffix}.csv" for suffix in SUFFIXES]
        for method in METHODS
    }
    canonical_names = [name for pair in names.values() for name in pair]
    prefix = "" if all((args.input_dir / name).is_file() for name in canonical_names) else "D_"
    paths = {name: args.input_dir / f"{prefix}{name}" for name in canonical_names}
    missing = [name for name, path in paths.items() if not path.is_file()]
    if missing:
        parser.error(f"{args.input_dir} must contain the complete Figure 6D CSV set "
                     f"(repository filenames or released D_ filenames): {', '.join(canonical_names)}")
    hashes = {name: hashlib.sha256(path.read_bytes()).hexdigest() for name, path in paths.items()}
    rows = turnover({method: [paths[name] for name in pair] for method, pair in names.items()})
    if any(hashlib.sha256(path.read_bytes()).hexdigest() != hashes[name] for name, path in paths.items()):
        raise RuntimeError("A Figure 6D input changed during calculation")

    args.output.mkdir(parents=True, exist_ok=False)
    table = args.output / "class_f1.csv"
    with table.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    manifest = {
        "paper_table": "Supplementary Table S11",
        "input_panel": "Figure 6D",
        "input_sha256": hashes,
        "output_sha256": {table.name: hashlib.sha256(table.read_bytes()).hexdigest()},
        "matched_genes": N_GENES,
        "gene_order": "lexicographic",
        "random_streams": "NumPy SeedSequence.spawn(2), bootstrap then permutation",
        "seed": SEED,
        "numpy": np.__version__,
        "bootstrap": {
            "draws": N_BOOTSTRAP,
            "strata": "scEU-seq reference class",
            "confidence_interval": "2.5th and 97.5th percentiles",
        },
        "permutation": {
            "draws": N_PERMUTATIONS,
            "scheme": "independent paired swaps of the two methods' labels within each gene",
            "statistic": "absolute class F1 difference",
            "batch_size": 2000,
            "tail_comparison_tolerance": 1e-12,
            "p_value": "(1 + number at least as extreme) / (1 + draws)",
            "adjustment": "Holm across the three classes",
        },
    }
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
