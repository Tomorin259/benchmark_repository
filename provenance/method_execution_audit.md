# Method execution audit

## Scope

This audit reconciles the 20 manuscript configurations with Supplementary Table S1, Figure 6, the 280-row cell-level render manifest, candidate method runners, final result paths, and surviving Conda-environment evidence. This pass reviewed source files, result paths, SHA-256 hashes, and import smoke tests; it did not rerun model training or regenerate scientific results.

The audit was completed on 2026-09-01 against the unversioned working root recorded in `phase1_baseline.md`, paper source `merged.tex` with SHA-256 `172bf2d84af52d8b7f005f588fbe81a44dc11f962bc4ccb5d6c35f7827b9f13a`, and release-repository base commit `7ad5c4d6afbf75707df93bd46d8a6d85720a0edc`.

## Confirmed configuration and result contracts

- The 16 top-level method descriptions expand to exactly 20 benchmark configurations.
- The render manifest and cell-level metric tables each contain 20 operational method identifiers across 14 datasets.
- After applying the 12 owner-confirmed TIvelo path aliases, every render-manifest result path resolves to an existing H5AD.
- Method-specific final velocity layers and embedding keys are recorded in `method_configurations.tsv`.
- The gene-level nine-method subset, Dynamo metabolic-labeling positive control, and RowGeometricMedian-4 utility are separate analysis roles and are not additional method configurations.
- For C10 (TFvelo), the official source is available at `/mnt/second19T/zhaozelin/TFvelo-main`. The chain is `TFvelo_run_demo.py` (preprocessing and recovery) followed by `TFvelo_analysis_demo.py` (analysis and final H5AD output). The new source-backed launcher records the 14 manual input overrides, stages the official relative paths, and patches only the analysis script's hard-coded hematopoiesis input in a temporary copy. The 14 files in `velocity:results/TFvelo/` are byte-identical to the corresponding files in the original TFvelo result tree.
- For C12 (veloVI), the cell-level canonical family is `velocity:results/veloVI/veloVI_without_vu/veloVI_*_all.h5ad` and stores the default spliced `get_velocity` output in `velocity`. The gene-level family is the root `velocity:results/veloVI/veloVI_*_all.h5ad` and additionally stores the explicit `velo_mode="unspliced"` output in `velocity_u`. The two roles are compatible with the observed code and results; the configuration status is confirmed.
- For C14 (TSvelo), the official source is available at `/mnt/second19T/zhaozelin/velocity/methods/TSvelo`. `TSvelo_run.py` calls the preprocessing and branch-model chain in `TSvelo_pp.py` before concatenating or renaming the final output. The source-backed launcher records all 14 manual input overrides and invokes the original entry point in an isolated compatibility tree. The mousebrain result intentionally uses the same biological dataset through a 4881-cell multiome loom, yielding 3653 cells after TSvelo preprocessing; cell-level evaluation aligns those cells to the common 3365-cell benchmark panel.

## Candidate runners, not yet canonical runners

The working project contains method-specific entry points for the 20 configurations. The registry records one runner per configuration and keeps source-backed launchers separate from the external upstream implementation. A manual method-local-to-final-result copy or rename is acceptable when its provenance is recorded for the released artifact. Historical path aliases and remaining environment limitations are recorded row by row in `runner_inventory.tsv`.

Known alternates excluded from the primary registry include partial or incompatible DeepKINET, cell2fate, GraphVelo, veloVI, veloAE, TSvelo, cellDancer, MultiVelo, and Dynamo scripts. They remain in the private working tree but are not silently promoted to public entry points.

For veloVI, the cell-level runner associated with the canonical `without_vu` result family calls `get_velocity(n_samples=25, velo_statistic="mean")`, whose velovi 0.3.1 default is `spliced`. The gene-level runner additionally calls `get_velocity(..., velo_mode="unspliced")` and stores `velocity_u`. The previously apparent manuscript/runner conflict is therefore a distinction between output roles, not an unresolved result identity.

## Final result update

The formal C13 result was subsequently frozen with benchmark seed 123. The
seeded runner resets Python, NumPy, and PyTorch/CUDA state before each dataset,
enables deterministic PyTorch algorithms, and preserves the upstream model
settings (`epoch=2000`, `loss_mode="poisson"`, `batch_key=None`). The 14
seed-123 H5AD outputs were reproduced in a shuffled dataset order with exact
core-layer equality; the only expected difference was the embedded checkpoint
path. The seed is a benchmark reproducibility control, not an upstream tutorial
default.

The C13 cell-level rows in the compact metric table and the RPE1/organoid
kinetics-recovery outputs were refreshed from this final result family. The
kinetics-recovery contract uses `degradation_rate` as DeepKINET gamma and
does not infer an alpha result that the H5AD does not provide.

## Environment evidence

The manuscript versions are well defined except for TFvelo, TSvelo, and DeepVelo_scv. Historical Conda assignments are mostly inferred. The import smoke audit found an importable candidate for all 16 configuration environments, but these checks did not execute model training. DeepVelo GCN imports with a NumPy/SciPy compatibility warning. The original SDEvelo environment fails during scVelo import because its Python 3.8/scVelo 0.3.3 combination is incompatible; the `sdevelo_py310` environment passes the smoke test and is now the C17 audited candidate. DeepKINET has a dirty nested repository and a specification that differs from the working environment; veloAE and DeepVelo_scv retain historical metadata or revision-cleanliness issues. The release therefore records runtime evidence separately from clean-spec readiness and does not generate misleading lock files.

## Remaining contract gaps

- Run the new TFvelo and TSvelo launchers on one ordinary and one special dataset, then compare their output contracts with the existing formal artifacts before any full rerun.
- Replace partial active lists with explicit supported-dataset mappings as each method is promoted; a shared universal runner is not required.
- Freeze the exact TSvelo preprocessing inputs, method revision, and bundled-asset licensing terms.
- Record manual curation and checksums for method-local outputs when those outputs are redistributed.
- Produce clean public environment specifications, including the audited SDEvelo candidate and the local-state cleanup required by DeepKINET.
- Resolve upstream revisions and redistribution terms before copying any method code or bundled databases.
- Rebuild the cell-metric input mapping because `all_metrics_main.py` currently has an empty `ADATA_PATHS` mapping.

## Release consequence

Phase 2 provides a reviewable execution registry, not a claim of all-20 one-command reproduction. The next code-curation phase should promote method-specific runners in small reviewed groups, beginning with entries whose dataset mappings and output contracts are explicit, while leaving remaining release work documented rather than guessed.
