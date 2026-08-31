# Method execution audit

## Scope

Phase 2 reconciles the 20 manuscript configurations with Supplementary Table S1, Figure 6, the 280-row cell-level render manifest, candidate method runners, final result paths, and surviving Conda-environment evidence. No scientific method was run and no large result or environment was copied.

The audit was completed on 2026-09-01 against the unversioned working root recorded in `phase1_baseline.md`, paper source `merged.tex` with SHA-256 `172bf2d84af52d8b7f005f588fbe81a44dc11f962bc4ccb5d6c35f7827b9f13a`, and release-repository base commit `ba9cf4abaaba1c3759831d6703653e29d85b6e26`.

## Confirmed configuration and result contracts

- The 16 top-level method descriptions expand to exactly 20 benchmark configurations.
- The render manifest and cell-level metric tables each contain 20 operational method identifiers across 14 datasets.
- After applying the 12 owner-confirmed TIvelo path aliases, every render-manifest result path resolves to an existing H5AD.
- Method-specific final velocity layers and embedding keys are recorded in `method_configurations.tsv`.
- The gene-level nine-method subset, Dynamo metabolic-labeling positive control, and RowGeometricMedian-4 utility are separate analysis roles and are not additional method configurations.

## Candidate runners, not yet canonical runners

The working project contains plausible entry points for 19 configurations, but many have partial active lists, dataset-specific wrappers, stale path aliases, or no frozen method-local-to-final-result copy step. The TFvelo result set exists while its local runner directory is empty. TSvelo has unresolved preprocessing provenance and duplicate result choices. These limitations are recorded row by row in `runner_inventory.tsv`.

Known alternates excluded from the primary registry include partial or incompatible DeepKINET, cell2fate, GraphVelo, veloVI, veloAE, TSvelo, cellDancer, MultiVelo, and Dynamo scripts. They remain in the private working tree but are not silently promoted to public entry points.

For veloVI, Supplementary Table S1 reports extraction with `velo_mode="unspliced"`, whereas the candidate runner associated with the canonical `without_vu` cell-level result family calls `get_velocity(n_samples=25, velo_statistic="mean")` without that argument. The audited velovi 0.3.1 signature defaults `velo_mode` to `spliced`. The result family remains confirmed as the benchmark input, but its exact generating invocation is therefore marked as conflicted pending reconciliation.

## Environment evidence

The manuscript versions are well defined except for TFvelo, TSvelo, and DeepVelo_scVelo. Historical Conda assignments are mostly inferred. No clean upstream commit is recorded for the Table S1 method set. DeepKINET has a dirty nested repository and conflicts between its checked specification and current environment; veloAE and DeepVelo_scVelo show conflicting installed metadata. The release therefore records environment evidence without generating misleading lock files.

## Remaining contract gaps

- Recover or reconstruct the TFvelo runner from an authoritative source.
- Reconcile the manuscript-reported veloVI velocity mode with the canonical cell-level result family and surviving runner evidence.
- Replace partial dataset lists with one manifest-driven invocation per configuration.
- Freeze the exact TSvelo preprocessing and duplicate-result selection.
- Record the curation link from method-local output to each final root result.
- Resolve upstream revisions and redistribution terms before copying any method code or bundled databases.
- Rebuild the cell-metric input mapping because `all_metrics_main.py` currently has an empty `ADATA_PATHS` mapping.

## Release consequence

Phase 2 provides a reviewable execution registry, not a claim of all-20 one-command reproduction. The next code-curation phase should promote configurations in small reviewed groups, beginning with runners whose 14-dataset mappings are already explicit, while leaving unresolved configurations documented rather than guessed.
