# RNA velocity benchmark reproduction repository

This repository is the curated release workspace for the RNA velocity benchmark. It is intentionally separate from the working project tree, which contains large datasets, method installations, checkpoints, exploratory analyses, and historical outputs.

## Current status

Phases 1 and 2 provide an audited scaffold. They record the paper-to-code/result mapping, data sources, the exact 20 method configurations, candidate execution entry points, historical environment evidence, and the intended public-release boundary. The repository is still not a one-command runnable release: Phase 2 records import smoke results separately from the work still required to prepare clean public specifications and runners.

No raw data, processed AnnData objects, trained checkpoints, or bulk method-result directories are included at this stage.

## Scientific scope

The paper evaluates 20 RNA velocity method configurations across:

- cell-level trajectory inference on 14 real datasets;
- gene-level velocity inference on simulated and real datasets;
- kinetic-parameter recovery on metabolic-labeling datasets;
- an optional post-hoc signed-transition ensemble, reported as a practical cell-level utility rather than a twenty-first RNA velocity model.

## Planned layout

```text
benchmark_repository/
├── preprocessing/                 # Dataset-specific preprocessing entry points
├── methods/                       # Curated method runners and adapters
├── envs/                          # Reproducible software environments
├── evaluation/
│   ├── cell/                      # Cell-level metrics
│   ├── gene/                      # Gene-level metrics
│   └── kinetic/                   # Kinetic-parameter metrics
├── simulation/                    # Simulation generation and evaluation
├── figure_reproduction/
│   ├── main_figures/
│   ├── supplementary_figures/
│   ├── plot_data/
│   └── manual_assets/
├── ensemble/                      # Post-hoc signed-transition ensemble
├── results_summary/               # Small, publication-facing result tables
├── manifests/                     # Source, release-scope, data, and result manifests
├── provenance/                    # Frozen source snapshots and audit notes
├── tests/                         # Smoke and numerical regression tests
└── benchmarking.md                # Detailed internal project navigation
```

## Phase 1 records

- `manifests/source_assets.tsv` maps every main and supplementary figure panel to the currently identified paper asset, generator or builder, direct source, and provenance status.
- `manifests/release_scope.tsv` states what will be included, curated, referenced externally, or excluded from the public repository.
- `provenance/phase1_baseline.md` freezes the source-tree and paper snapshot used for this inventory.
- `DECISIONS.md` lists unresolved publication decisions without silently choosing them.
- [`manifests/datasets.tsv`](manifests/datasets.tsv), [`manifests/publication_artifacts.tsv`](manifests/publication_artifacts.tsv), and [`manifests/path_aliases.tsv`](manifests/path_aliases.tsv) record the audited source and artifact contract.
- [`provenance/data_source_audit.md`](provenance/data_source_audit.md) summarizes evidence and release consequences.

## Phase 2 records

- [`manifests/method_configurations.tsv`](manifests/method_configurations.tsv) freezes the manuscript-to-operational identity of the 20 configurations and their final output keys.
- [`manifests/runner_inventory.tsv`](manifests/runner_inventory.tsv) records candidate runners, coverage limitations, and observed result contracts without promoting partial scripts to canonical entry points. A method may have its own release script; a universal runner is not required.
- [`manifests/environments.tsv`](manifests/environments.tsv) separates manuscript-reported versions and provenance evidence from `runtime_status` import smoke results and `release_spec_status` public-specification readiness.
- [`methods/README.md`](methods/README.md), [`envs/README.md`](envs/README.md), and [`provenance/method_execution_audit.md`](provenance/method_execution_audit.md) define the code-curation boundary.
- `python tests/validate_phase2_manifests.py` checks the registry contract using only the Python standard library.

## Provenance policy

The working project root is not version-controlled as a whole. Therefore, a path or modification time alone is not treated as proof of provenance. A publication-facing artifact must be tied to at least one of the following:

1. a Git commit;
2. a file hash and frozen manifest;
3. a direct result manifest used by the plotting code;
4. a documented manual assembly step backed by canonical source assets.

Historical folders are excluded by default, except when a current paper panel demonstrably uses an asset stored there (currently Figure 5A). Manual result copying or renaming is recorded as provenance when relevant; it is not treated as an automatic-runner requirement.

## Not yet provided

This Phase 2 registry is not a runnable release. Reviewed environment creation, promotion of thin method runners, metric reproduction, figure regeneration, data-access instructions, license files, and archival deposition are later review gates.
