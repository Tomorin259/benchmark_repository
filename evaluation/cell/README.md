# Cell-level evaluation contract

`config.json` is the Phase 3A machine-readable extraction of the
current effective `indicators/cell/all_metrics_main.py` source.  The source
symbols used for the extraction are recorded in the JSON metadata; this file
does not copy or depend on a private absolute path.

## Dataset configuration

The `datasets` object has exactly 14 keys, preserving the dataset names used by
`DATASET_CONFIGS` and by the current metric CSV:

| source key | `dataset_id` foreign key | `emb` | `label` | BATC `n_bins` |
| --- | --- | --- | --- | ---: |
| `bonemarrow` | `bonemarrow` | `umap` | `clusters` | 30 |
| `capofibro` | `capo_pre` | `umap` | `clusters` | 5 |
| `DentateGyrus_lamanno` | `dentategyrus_lamanno` | `umap` | `clusters` | 20 |
| `dentategyrus` | `dentategyrus` | `umap` | `clusters` | 30 |
| `forebrain` | `forebrain` | `umap` | `clusters` | 20 |
| `gastrulation` | `gastrulation` | `umap` | `clusters` | 20 |
| `hematopoiesis` | `hematopoiesis` | `umap` | `clusters` | 20 |
| `mousebrain` | `mousebrain` | `umap` | `clusters` | 20 |
| `neuron` | `scnt_neuron` | `umap` | `clusters` | 20 |
| `organoid` | `sceu_organoid` | `umap` | `cell_type` | 40 |
| `pancreas` | `pancreas` | `umap` | `clusters` | 30 |
| `retina` | `retina` | `umap` | `clusters` | 40 |
| `rpe1` | `sceu_rpe1` | `umap` | `clusters` | 10 |
| `rpe1_gen` | `rpe1_gen` | `umap` | `clusters` | 20 |

The `dataset_id` values must resolve to the cell-level rows in
`data/datasets.tsv`.  `rpe1_gen` is not present in the source
`BATC_N_BINS_BY_DATASET` map, so its effective value is the extracted
`DEFAULT_BATC_N_BINS` value of 20.  Every `cluster_edges` entry is a two-item
array; string and integer endpoints are preserved exactly, including the
integer endpoints for `forebrain`.

### Label-source behavior

For the ordinary `clusters` datasets, the generic fallback order is exactly:

```text
celltype -> Clusters -> phase -> Annotation -> cell_type
```

It is consulted only when `clusters` is absent.  Before that fallback,
`prepare_clusters` copies the following dataset-specific column when it is
present and then returns:

| dataset | source column | target |
| --- | --- | --- |
| `rpe1` | `cell_cycle_phase` | `clusters` |
| `retina` | `Annotation` | `clusters` |
| `neuron` | `time_str` | `clusters` |
| `hematopoiesis` | `cell_type` | `clusters` |
| `forebrain` | `Clusters` | `clusters` |
| `gastrulation` | `celltype` | `clusters` |

For `organoid`, `prepare_clusters` returns immediately and the configured
`cell_type` column must already exist.  `CELL_MAPPER_ORGANOID` and its
`som_cluster_id` mapping are recorded for provenance with the literal status
`legacy_defined_but_not_applied`; the mapper is not an executed preparation
step and must not be described as one.

## Method embedding overrides

The static `METHOD_EMB_OVERRIDES` extraction contains all 20 operational
methods.  The source normalizes method names to lowercase before lookup, so
the validator compares the JSON keys case-insensitively.  All methods use
`umap` except `dynamo`, which uses `umap_ori`:

```text
cell2fate, celldancer, deepkinet, deepvelo, deepvelo2022,
graphvelo_dynamo, graphvelo_scv, multivelo, scvelo_deterministic,
scvelo_dynamic, scvelo_stochastic, sdevelo, TFvelo, tivelo, tsvelo,
unitvelo, unitvelo_mode2, veloae, velovi -> umap
dynamo -> umap_ori
```

This is a method-level override.  An explicit command-line embedding override
in the source remains a higher-priority runtime input and is not represented
as a second default in this contract.

## Frozen result-input manifest

The private audit file `manifests/cell_result_inputs.tsv` uses this exact header
and contains one row
for every `(config_id, dataset_id)` pair (20 × 14 = 280 rows):

```text
config_id\toperational_id\tdataset_id\tevaluation_dataset\tsource_result_path\tresult_path\tpath_resolution\tvelocity_layer_key\tsource_embedding_key\tvelocity_embedding_key\tevaluation_embedding\talignment_mode\tn_reference\tn_method\tn_common\tsource_embedding_max_error\taffine_relative_rmse\tsource_nonfinite_velocity_values\taligned_nonfinite_velocity_values\tfile_size_bytes\tsha256
```

`config_id` and `operational_id` are foreign keys to
`methods/configurations.tsv`; `dataset_id` is a foreign key to
`data/datasets.tsv`. `evaluation_dataset` and `operational_id` form the
key used by the copied metric table. The three layer/embedding fields and the
effective evaluation embedding were inspected against the original H5AD files
during the Phase 3A audit.

`source_result_path` preserves the path stored by the source render manifest.
`result_path` is a private audit locator. The public release replaces it with
an artifact ID from the root `artifacts.tsv`; absolute or server-specific paths
are not part of the public contract. `file_size_bytes` is a positive decimal
integer and `sha256` is the digest of the canonical file bytes. Alignment fields
are evidence copied from the render manifest, not newly calculated trajectory
metrics.

## Metric table

`results/all_metrics.csv` uses the source
header:

```text
dataset,method,consistency,cbdir,iccoh,batc
```

Its 280 unique `(dataset, method)` keys must equal
`(evaluation_dataset, operational_id)` from `cell_result_inputs.tsv`.
`consistency`, `cbdir`, `iccoh`, and `batc` must all
parse as finite floating-point values. The file is copied byte for byte from
the audited working result and was not recalculated in Phase 3A.

Run the standard-library validator from the repository root:

```bash
python tests/validate_cell_contract.py
```

Use `--require-artifacts` for the release/CI gate. To check paths and recorded
file sizes against a local working tree, add `--working-root PATH`; add
`--verify-hashes` to compute and compare SHA-256 as well.  No H5AD or AnnData
inspection is performed by the standard-library validator; the backed,
read-only key audit is recorded separately in the provenance report.
