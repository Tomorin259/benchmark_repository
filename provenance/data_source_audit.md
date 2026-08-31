# Data-source audit

## Scope and source of record

This audit merges `merged.tex`, the remote working tree, and the scVelo 0.3.4 source snapshot.

## Confirmed

The 17-dataset inventory is complete. Exact upstream copies are confirmed for Pancreas, Gastrulation, Dentategyrus, and Bonemarrow; the source loom is exact for Mousebrain; Hematopoiesis v2 matches by MD5. Canonical result/table chains are identified.

## Inferred/derived

Several H5ADs are processed subsets or transformations. Simulation inputs are generated artifacts, and final-gene collections are manually assembled evidence.

## Reproducibility gaps

Retina, Mousebrain H5AD, scNT, capo_pre, and rpe1_gen transforms are not frozen. Forebrain/DG legacy endpoints return HTML and cannot be byte-checked. Simulation seed/output path/revision metadata is missing. Cell metric CSVs have no per-row source and `ADATA_PATHS` is empty. DG_lamanno has two H5AD variants. Some method runners are stale; do not claim identical filenames.

## Release consequence

Git gets scripts, manifests, and compact tables; canonical processed H5ADs and selected outputs go to a checksum-addressed external archive. Do not claim raw-to-result one-command reproduction.

## Verification table

| Input | Verification |
|---|---|
| Pancreas | SHA match |
| Gastrulation | SHA match |
| Dentategyrus | SHA match |
| Bonemarrow | SHA match |
| Hematopoiesis v2 | MD5 exact |
| Mousebrain source loom | SHA exact |

