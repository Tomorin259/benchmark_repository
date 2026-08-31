# Data-source audit

## Scope and source of record

This Phase 1 audit reconciles the 17 datasets named in `merged.tex` with the working project at `velocity:`, `library.bib`, scVelo 0.3.4 source paths, publication result manifests, and authoritative upstream records. It distinguishes an official biological source from the exact processed file consumed by the benchmark: an accession alone does not prove the local preprocessing chain.

The inventory contains 14 real datasets and three generated simulations. Audited dimensions and SHA-256 values describe the local benchmark inputs, not necessarily the dimensions of the upstream raw dataset.

## Exact source-file verification

Six upstream relations were verified at byte or checksum level:

| Audited input | Evidence |
|---|---|
| Pancreas | Local SHA-256 equals the scVelo tutorial H5AD. |
| Gastrulation | Local SHA-256 equals Figshare file `27686871`. |
| Dentategyrus | Local SHA-256 equals the scVelo tutorial H5AD. |
| Bonemarrow | Local SHA-256 equals Figshare file `27686835`. |
| Hematopoiesis | Local MD5 `c738b0ea2246d733440ca11b2b87674d` equals Dynamo Figshare v2 file `47439626`. |
| Mousebrain source loom | Local SHA-256 equals public file `54153947`; the benchmark H5AD derived from it is not byte-identical and its transform is not frozen. |

All exact hashes and URLs are recorded in `manifests/datasets.tsv`.

## Confirmed source with unresolved local transformation

- **Forebrain:** SRP129388 and the scVelo/velocyto loom endpoint are confirmed, but the endpoint now returns HTML rather than the historical loom. The local loom must therefore be frozen by checksum in the external archive.
- **Dentategyrus_lamanno:** GSE104323 is the developing-dentate-gyrus source used by La Manno et al. The historical loom endpoint and conversion code are identifiable, but endpoint byte identity is no longer testable. Both benchmark H5AD variants remain in scope because different runners reference them.
- **Retina:** GSM3466902/GSE122466 is the source record. The audited input has 2,726 cells, whereas the GEO sample title reports 2,675; the exact local transformation is missing.
- **Mousebrain:** the 10x E18 source loom is exact, but the 3,365-cell, 936-gene benchmark H5AD transformation is not frozen.
- **scEU_rpe1 and scEU_organoid:** the biological sources and corresponding Dynamo Figshare v2 processed files are identified; local serialization or downstream transformation is not fully reconstructed.
- **capo_pre and rpe1_gen:** the VeloCycle Zenodo files have the same cell identifiers as the local inputs. The local files are 2,000-gene processed subsets, and the exact subset/preprocessing commands are missing.
- **scNT_neuron:** the Dynamo neuron-splicing source is the best-supported origin, but the local 3,066-cell excitatory time-course subset and 2,000-gene transform remain an inference.

These inputs are valid frozen benchmark artifacts, but the release must not describe them as one-command raw-to-input reproductions until the missing transforms are recovered.

## Simulations

The canonical generated inputs and their structurally corresponding generators are listed in `manifests/datasets.tsv`. Array-level checks verified the stored truth fields against the generator equations with maximum absolute error zero:

\[
v_u = \alpha - \beta u, \qquad v_s = \beta u - \gamma s.
\]

The scripts do not fix the simulation random seed, their output paths are stale, and no original invocation or source revision is embedded in the H5ADs. The canonical files therefore require external archival checksums; the present scripts establish the model structure, not byte-exact regeneration.

## Publication result chains

- `velocity:plot/unified_velocity_stream/results/render_manifest.csv` contains the 14-dataset by 20-method stream lookup. Twelve TIvelo entries retain historical paths; `manifests/path_aliases.tsv` records the owner-confirmed moves to existing root-level TIvelo files.
- `all_metrics.csv` and `transprob.csv` are the compact canonical cell-level tables. Their present rows do not retain the H5AD input path, and the current `all_metrics_main.py` has an empty `ADATA_PATHS` mapping. The tables can be published, but exact row-level recomputation needs a new frozen input manifest.
- Current simulation scores, real-gene score tables, Figure 6 summary CSVs, and the RowGeometricMedian-4 compact result package are identified in `manifests/publication_artifacts.tsv`.
- Some final figures were manually composed from verified current collections. This is documented as a packaging gap, not silently treated as an automated figure build.

## Release consequence

Phase 1 determines the source and canonical local artifact for every dataset; it does not claim that every transformation is reproducible. The public Git repository should contain scripts, manifests, compact numerical tables, and selected small figure inputs. Canonical processed H5ADs, simulation inputs, and any required large method outputs should be deposited as a checksum-addressed external archive. The unresolved transformation, licensing, and archival-host decisions remain explicit review gates in `DECISIONS.md`.
