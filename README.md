# RNA velocity benchmark reproduction repository

This repository curates the code and compact numerical outputs needed to
understand and reproduce the RNA velocity benchmark. Large datasets, method
result H5AD files, checkpoints, environments, logs, and exploratory outputs
remain in external storage.

## Current status

As of 2026-09-04, the paper-aligned repository skeleton is complete. The only
evaluation module currently implemented as runnable code is
`evaluation/parameter_recovery/`. The cell, gene, summary, method-reproduction,
environment, and ensemble areas currently contain only the audited material or
structural placeholders needed for their later curation; their presence must
not be read as implementation completion.

## Repository path

```text
methods/                         method-specific reproduction scripts
envs/                            reviewed software specifications
evaluation/cell/                 cell-level evaluation and panels
evaluation/gene/                 gene-level evaluation and simulation
evaluation/parameter_recovery/  kinetic-parameter recovery (paper Case 3)
evaluation/summary/              cross-experiment summaries such as Figure 6
ensemble/                        optional post-hoc signed-transition utility
data/                            dataset metadata and access information
artifacts.tsv                    portable index of external artifacts
docs/                            cross-module documentation and figure map
tests/                           lightweight publication-contract checks
```

There is intentionally no shared top-level preprocessing package, results
directory, or figures directory. Method-specific preprocessing stays in the
method script that requires it; evaluation processing, compact results, and
generated panels stay together in the evaluation module that owns the
scientific question.

## Running parameter recovery

The parameter-recovery module consumes saved method-result H5AD files and the
scEU-seq reference tables. It does not train methods or prepare the benchmark
datasets.

```bash
python evaluation/parameter_recovery/metrics.py \
  --artifact-root /path/to/artifacts \
  --check-only

python evaluation/parameter_recovery/metrics.py \
  --artifact-root /path/to/artifacts \
  --output /path/to/metrics

python evaluation/parameter_recovery/plots.py \
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
cell2fate will document the four special datasets, cell/gene alignment, finite
and non-negative count checks, integer-valued checks, and lossless dtype
conversion inside their own method scripts when those scripts are curated.
Methods that do not require counts do not use this path.

## Private audit material

The current private repository also retains `manifests/` and `provenance/` as
the evidence trail for curation. These directories are not part of the clean
public release interface. The public-facing artifact index is the root
`artifacts.tsv`; it must not expose server-specific working paths.

Read [`AGENTS.md`](AGENTS.md) and [`REPOSITORY_DESIGN.md`](REPOSITORY_DESIGN.md)
before changing a script or module.
