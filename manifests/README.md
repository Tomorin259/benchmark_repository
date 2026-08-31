# Manifest conventions

All paths in the Phase 1 manifests are relative to one of two roots:

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

`datasets.tsv` is the dataset contract: identity, dimensions, references, source URLs, audited local inputs, checksums, relation, evidence status, and release action. `publication_artifacts.tsv` maps publication consumers to source paths and storage policy. `path_aliases.tsv` records historical-to-current path drift and evidence. All are tab-separated and use `velocity:` locators rather than absolute paths. Official source URLs identify authoritative records; audited local inputs are files inspected in the working tree. `confirmed`, `unresolved`, `mixed`, and `inference` distinguish byte-level evidence, gaps, combined source/derived content, and structural evidence.
