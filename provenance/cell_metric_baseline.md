# Cell-level evaluation baseline

## Scope

Phase 3A freezes the inputs and existing outputs of the cell-level benchmark. It does not rerun any RNA velocity method, recompute a metric, or change the evaluation algorithm. The audit was performed on 2026-09-01 from release-repository commit `9f56bd6b11ed1eb28c9bca29ab95ed27386893e0`.

The contract grain is one row per method configuration and evaluation dataset: 20 configurations by 14 datasets, for 280 unique `(config_id, dataset_id)` rows. The copied metric table has the corresponding 280 unique `(dataset, method)` rows used by the working evaluation output.

## Frozen sources

| Artifact | Audited source | SHA-256 |
| --- | --- | --- |
| Render manifest | `velocity:plot/unified_velocity_stream/results/render_manifest.csv` | `f66a636e9629b25c110a18d96c28280fcb3b4c3566458e077641b28bb471eecf` |
| Metric baseline | `velocity:indicators/cell/all_metrics/all_metrics.csv` | `f08e6f3f32c179541b80dda24c7cb2c7a0434701c929fa92f02c8073cb97faff` |
| Legacy metric source | `velocity:indicators/cell/all_metrics_main.py` | `c79929214882afc5d8dcb5e219b8a8f16785eb5a7938143f6e1e5c43d1885c50` |
| Frozen input manifest | `manifests/cell_result_inputs.tsv` | `bd639c5f34f7cc4fefedd093f927babf43b12d3f754b69daf1638be417666f3d` |

`results_summary/cell_metrics/all_metrics.csv` is a byte-identical copy of the audited metric baseline: its SHA-256 is also `f08e6f3f32c179541b80dda24c7cb2c7a0434701c929fa92f02c8073cb97faff`. Negative CBDir and BATC values are valid signed scores and were not filtered or clipped.

## Canonical result resolution

The source render manifest contains 12 historical TIvelo paths under `TIvelo_for_cells`. Those files had been moved before this audit. For each affected row, `source_result_path` preserves the render-manifest path, while `result_path` stores the owner-confirmed current path from `manifests/path_aliases.tsv` and `path_resolution` is `owner_confirmed_move`. The other 268 rows are `direct`.

All 280 canonical paths exist and are unique. Their combined size is 227,223,730,170 bytes (211.62 GiB). Each canonical file was read once with GNU `sha256sum`, using eight file-level workers; the resulting 280 SHA-256 digests are all unique and are stored with the exact file sizes in `cell_result_inputs.tsv`.

## Structural and H5AD checks

The following checks passed:

- exactly 20 configurations and 14 datasets, with 14 rows per configuration and 20 rows per dataset;
- no duplicate `(config_id, dataset_id)` key and no duplicate canonical result path;
- the metric-table key set exactly matches `(evaluation_dataset, operational_id)` from the input manifest;
- all four metric columns parse as finite floating-point values;
- all 280 recorded file sizes match the canonical files;
- all 280 H5AD files open read-only;
- every recorded velocity layer, source embedding and velocity embedding key exists in its H5AD file.

The H5AD metadata audit used Python 3.10.19 and h5py 3.16.0 from the audited `velocity_ensemble` environment. It opened only HDF5 metadata and did not load matrices or modify files.

## Preserved alignment evidence

The render-manifest evidence is retained rather than normalized away:

- alignment modes: 271 `direct`, 6 `graph_reproject`, and 3 `affine` rows;
- five rows have method/reference cell-count differences:
  - TFvelo–bonemarrow: 5,780 reference, 5,756 method/common cells;
  - Dynamo–dentategyrus: 2,930 reference, 2,918 method/common cells;
  - TFvelo–hematopoiesis: 1,947 reference, 1,773 method/common cells;
  - TSvelo–mousebrain: 3,365 reference, 3,653 method cells and 3,365 common cells, using graph reprojection;
  - Dynamo–retina: 2,726 reference, 2,683 method/common cells;
- TSvelo–bonemarrow has 986 non-finite values in the source velocity embedding and zero after alignment; no other row has a nonzero source or aligned non-finite count.

These observations are provenance for the existing plotted/evaluated results. They are not treated as missing rows and do not establish that the underlying method runs have been independently reproduced.

## Dataset configuration boundary

`evaluation/cell/dataset_config.json` is a machine-readable extraction of the effective dataset priors and routing in the legacy metric source: label columns, directed cluster edges, BATC bin counts, default embeddings, label-source rules and method embedding overrides.

The legacy source defines `CELL_MAPPER_ORGANOID`, but `prepare_clusters` returns immediately for organoid and does not apply that mapping. The JSON therefore records it as `legacy_defined_but_not_applied`; the existing `cell_type` column is the effective input. This distinction must be preserved when the metric code is curated in Phase 3B.

## Limits

This freeze proves the identity, tabular completeness, file checksums, declared H5AD keys and copied metric values of the current cell-level result set. It does not prove that the 280 H5AD files can yet be regenerated from clean method environments, that the historical metric values have been independently recomputed, or that all stochastic methods will reproduce byte-identical outputs. Those are later execution and numerical-regression gates.
