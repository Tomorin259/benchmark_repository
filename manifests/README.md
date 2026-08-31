# Manifest conventions

All paths in the Phase 1 and Phase 2 manifests are relative to one of two roots:

- `velocity:` — `/mnt/second19T/zhaozelin/velocity` in the audited working tree;
- `paper:` — the frozen paper project snapshot identified in `provenance/phase1_baseline.md`.

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

`method_configurations.tsv` is the exact 20-row manuscript registry. `operational_id` follows the render and metric tables; `manuscript_label` preserves Supplementary Table S1 naming. A `confirmed` configuration status confirms its identity and observed final-output contract, not that its runner is publication-ready. `conflicted` retains a confirmed result identity while flagging a material mismatch between the manuscript setting and surviving execution evidence.

`runner_inventory.tsv` records one audited candidate row per configuration. Status values are `candidate`, `partial`, `ambiguous`, and `missing`. Candidate and partial scripts require review before copying; `recover_before_release` means the execution chain must be recovered rather than inferred from outputs.

`environments.tsv` records one row per distinct audited Conda assignment. Absolute paths are allowed only in `audited_conda_path` because they are private provenance evidence, not public commands. Status values distinguish `partial`, `ambiguous`, `conflicted`, and `missing` evidence. These rows must not be converted mechanically into lock files.

Run `python tests/validate_phase2_manifests.py` after editing any Phase 2 registry.
