# RNA velocity benchmark reproduction repository

This repository is the curated release workspace for the RNA velocity benchmark. It is intentionally separate from the working project tree, which contains large datasets, method installations, checkpoints, exploratory analyses, and historical outputs.

## Current status

Phase 1 is an inventory-only scaffold. It records the paper-to-code/result mapping, data-source audit, and intended public-release boundary. It is still not runnable.

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

## Provenance policy

The working project root is not version-controlled as a whole. Therefore, a path or modification time alone is not treated as proof of provenance. A publication-facing artifact must be tied to at least one of the following:

1. a Git commit;
2. a file hash and frozen manifest;
3. a direct result manifest used by the plotting code;
4. a documented manual assembly step backed by canonical source assets.

Historical folders are excluded by default, except when a current paper panel demonstrably uses an asset stored there (currently Figure 5A).

## Not yet provided

This Phase 1 commit is not a runnable release. Environment creation, method execution, metric reproduction, figure regeneration, data-access instructions, license files, and archival deposition are later review gates.
