# RNA velocity benchmark reproduction repository

This repository curates the code and compact numerical outputs needed to
understand and reproduce the RNA velocity benchmark. The benchmark has 20
method configurations, registered only in
methods/configurations.tsv. Large datasets, method result H5AD files,
checkpoints, environments, logs, and exploratory outputs remain in external
storage.

The two post-hoc comparators are separate from those 20 methods: Signed-T
Mean-4 is a baseline and Signed-T RowGeometricMedian-4 is the formal ensemble.
Each is ranked against
the 20 methods separately; the two comparators are not ranked against each
other. Mean-4 is the arithmetic mean of the four signed transition values on
each common edge. RowGeometricMedian-4 takes the Euclidean geometric median of
the four complete signed off-diagonal transition rows independently for each
source cell, without post-hoc L1 renormalization. For each dataset and metric,
each variant is ranked independently
against the 20 Cell methods: rank is 1 plus the number of strictly higher method
scores, with exact method ties reported using competition-rank tie counts; the
ensemble is participant 21 and variants are never candidates for one another.
Their implementations are in ensemble/, using the portable fixed-neighbor
assets in evaluation/cell/reference_neighbors/. Those assets install the
audited candidate graph and indices, not their build-time PCA coordinates;
cosines are still computed in each input method's native Ms/velocity space.
Formal fixed-neighbor all-14 ensemble results are
complete and registered in `artifacts.tsv`. The canonical scores, summaries, and
independent rankings are [`ensemble/results/scores.csv`](ensemble/results/scores.csv),
[`ensemble/results/summary.md`](ensemble/results/summary.md),
[`ensemble/results/rankings_vs_20.csv`](ensemble/results/rankings_vs_20.csv), and
[`ensemble/results/ranking_summary.md`](ensemble/results/ranking_summary.md).

## Current status

As of 2026-09-05, the paper-aligned repository contains the final seed-123
DeepKINET runner, a runnable cell-level evaluator, refreshed compact tables,
and a runnable kinetics-recovery bundle for the current DeepKINET result. The
gene and integrated summary evaluators remain structural placeholders or
curated data packages. The presence of a directory must not be read as full
all-20 benchmark completion.

## Repository path

```text
methods/                         method-specific reproduction scripts
envs/                            reviewed software specifications
evaluation/cell/                 cell-level evaluation and panels
evaluation/gene/                 gene-level evaluation and simulation
evaluation/kinetics_recovery/  kinetic-parameter recovery and its plotting module
evaluation/summary/              cross-experiment summaries such as Figure 6
ensemble/                        optional post-hoc signed-transition utility
data/                            dataset metadata and access information
artifacts.tsv                    portable index of external artifacts
figures/                         confirmed paper-level composite PDFs
docs/                            cross-module documentation and figure map
tests/                           lightweight publication-contract checks
```

There is intentionally no shared top-level preprocessing package or results
directory. The top-level figures directory contains only confirmed paper-level
composite PDFs; panel assets stay in the evaluation module that owns the
scientific question.

Method names shown in public tables use veloVI spelling. Dataset IDs identify
the biological input and evaluation panel, not methods or ensemble variants.

## Running kinetics-recovery

The kinetics-recovery module consumes saved method-result H5AD files and the
scEU-seq reference tables. It does not train methods or prepare the benchmark
datasets.

```bash
python evaluation/kinetics_recovery/metrics.py \
  --artifact-root /path/to/artifacts \
  --check-only

python evaluation/kinetics_recovery/metrics.py \
  --artifact-root /path/to/artifacts \
  --output /path/to/metrics

python evaluation/kinetics_recovery/plotting/figure5.py \
  --metrics /path/to/metrics \
  --output /path/to/figures
```

The evaluator extracts each method's native parameter keys at its input
boundary, applies the fixed RPE1/organoid processing rules, calculates
per-gene and summary metrics, and writes compact CSV files. Plotting consumes
only that metrics bundle.

## Counts-specific inputs

Only methods that genuinely require raw counts select the corresponding
`_counts.h5ad` artifact and perform its method-specific checks. DeepKINET and
cell2fate document the four special datasets, cell/gene alignment, finite and
non-negative count checks, integer-valued checks, and lossless dtype conversion
inside their own method scripts. DeepKINET additionally records the seed-123
deterministic benchmark setting used by the formal result.
Methods that do not require counts do not use this path.

## Public artifact index

The public-facing artifact index is the root artifacts.tsv; it must not expose server-specific working paths.

Read AGENTS.md and REPOSITORY_DESIGN.md before changing a script or module.
