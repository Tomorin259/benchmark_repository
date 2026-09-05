# Manifest conventions

All paths in the Phase 1 and Phase 2 manifests are relative to one of two roots:

- `velocity:` — `/mnt/second19T/zhaozelin/velocity` in the audited working tree;
- `paper:` — the frozen paper project snapshot identified in `provenance/phase1_baseline.md`.

The `project:` locator is used for audited source trees that sit beside the working
project: `project:` maps to `/mnt/second19T/zhaozelin` in the private audit
environment. It is a provenance locator only and is not a public filesystem path.

These prefixes are provenance locators, not intended public paths. Later phases will replace them with repository-relative paths or external accession URLs.

## `source_assets.tsv`

Each row describes one paper panel or one unsplit figure. Repeated paper hashes are expected because several panels are embedded in one final PDF.

Status values:

- `confirmed_current`: the named direct source and generator/builder are present and correspond to the current paper chain;
- `confirmed_manual_assembly`: canonical direct sources are present, but the final PDF was assembled or cropped manually;
- `confirmed_collection_manual_assembly`: the canonical result collection is identified, while exact panel packaging remains a later figure-release task;
- `manual_source_missing`: the final paper asset exists, but no editable source or generator was found;
- `builder_missing`: direct numerical inputs are known, but no complete final-figure builder was found.

`confirmed` means supported by the current file tree, paper source, a result manifest, numerical comparison, or Git history. It does not imply that the full workflow has already been rerun in the new release repository.

## `release_scope.tsv`

Decision values:

- `include`: copy a small, canonical artifact as-is after verification;
- `curate`: select and adapt only the publication-relevant subset;
- `external`: provide accession/checksum or upstream source instead of copying large/bundled content;
- `exclude`: do not publish the working-tree material;
- `decision_needed`: wait for an explicit scientific or licensing decision.

## Audit schemas

`datasets.tsv` is the dataset contract: identity, audited local dimensions, references, source URLs, audited local inputs, checksums, relation, evidence status, and release action. `publication_artifacts.tsv` maps publication consumers to source paths and storage policy. `path_aliases.tsv` records historical-to-current path drift and evidence. All are tab-separated and use `velocity:` locators rather than absolute paths. Official source URLs identify authoritative records; audited local inputs are files inspected in the working tree. `confirmed` denotes an exact or otherwise fully established relation, `mixed` denotes a confirmed source with an unresolved local transformation, `inference` denotes a structurally supported but unproven relation, and `unresolved` is reserved for a material gap without sufficient evidence.

For `publication_artifacts.tsv`, `git` means a compact artifact should be copied into the release repository, `curated` means only the named publication-facing subset should be copied, `external_archive` means the artifact should be frozen outside Git with checksums, and `decision_needed` marks a licensing or release-boundary question that remains open.

## Phase 2 execution schemas

`method_configurations.tsv` is the exact 20-row manuscript registry. `operational_id` follows the render and metric tables; `manuscript_label` preserves Supplementary Table S1 naming. A `confirmed` configuration status confirms its identity and observed final-output contract, not that its runner is publication-ready. The configuration row may describe more than one output role when the method has separate cell-level and gene-level outputs; veloVI records its spliced cell-level output and its additional unspliced gene-level output in this way.

`runner_inventory.tsv` records one audited candidate row per configuration. Status values are `candidate`, `source_backed`, `partial`, `ambiguous`, and `missing`. `source_backed` identifies a repository launcher that invokes an external official source checkout while recording its input contract. `missing` is reserved for a configuration for which no authoritative execution source has been located; an empty method subdirectory does not imply missing source when an external project locator is documented. Candidate, source-backed, and partial scripts require curation before publication. `recover_before_release` is reserved for a genuinely absent execution chain; a partial method-specific script uses `curate_after_review`. A final-result copy or rename may be manual and is recorded as provenance when the artifact is released.

`environments.tsv` records one row per distinct audited Conda assignment. Absolute paths are allowed only in `audited_conda_path` because they are private provenance evidence, not public commands. `evidence_status` describes the strength or limitations of the historical assignment evidence. `runtime_status` records only the observed package/entry-point import smoke result, with values `import_smoke_passed`, `import_smoke_passed_with_warning`, and `import_smoke_failed`; it is not evidence of a completed model run. `release_spec_status` records whether a clean public environment specification has been curated, with values `not_yet_curated`, `curation_required`, and `reviewed_clean`. These fields are intentionally independent: an environment can pass import smoke while still requiring release curation. These rows must not be converted mechanically into lock files.

Run `python tests/validate_phase2_manifests.py` after editing any Phase 2 registry.

## Phase 3A cell-level result manifest

The publication-facing cell evaluator consumes
`evaluation/cell/result_inputs.tsv`, whose paths are relative to an external
artifact root. The table below describes the separate private audit manifest;
it is retained for result provenance and contract checks, not as the runtime
evaluation input.

`manifests/cell_result_inputs.tsv` is the frozen 280-row (20 configurations ×
14 datasets) result-input manifest. Its primary key is `(config_id,
dataset_id)`, with `operational_id` and `evaluation_dataset` retaining the
names used by the published metric table. `source_result_path` preserves the
render-manifest path, while `result_path` identifies the current canonical
file. Twelve owner-confirmed TIvelo moves are represented explicitly by
`path_resolution=owner_confirmed_move`; all other rows are `direct`.

Each row also records the velocity layer, source and velocity embedding keys,
effective evaluation embedding, alignment mode and dimensions, non-finite
value diagnostics, positive file size, and a 64-character SHA-256 digest. The
exact schema and allowed `velocity:` locators are documented in
[`evaluation/cell/README.md`](../evaluation/cell/README.md).
The companion `evaluation/cell/results/all_metrics.csv` keeps the source
header `dataset,method,consistency,cbdir,iccoh,batc` and has the same 280-key
set with four finite metrics. `evaluation/cell/results/transprob.csv` is the
companion 280-row TransProb table. These are compact publication results,
whereas the audit manifest remains a private locator record.

Run `python manifests/validate_cell_contract.py` for the standard-library checks.
It does not open H5AD/AnnData files.  Pass `--working-root PATH` for optional
file existence and size checks, and add `--verify-hashes` for SHA-256 checks;
the default invocation has no private working-tree path dependency.

## Kinetics-recovery

evaluation/kinetics_recovery/config.json is the direct input map for the
kinetic-parameter recovery experiment. It contains reference-table paths,
canonical method H5AD paths, parameter-layer mappings, and the explicit
Dynamo/TSvelo H5AD extraction rules. Method
execution and dataset preparation remain outside this release layer; historical
method-specific CSVs are regression baselines rather than evaluator inputs.

evaluation/kinetics_recovery/metrics.py is the single metrics entry point. It writes
processed trajectories, per-gene metrics, component means, and equal-system
summaries to a new output directory. Use --check-only to validate direct
inputs before producing metrics. Future Figure 5/6 scripts should consume
that output instead of reopening unlisted result files.
