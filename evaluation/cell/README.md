# Cell-level evaluation

This module evaluates the 20 configured RNA velocity methods across the 14
real benchmark datasets. It consumes method-native H5AD files from an external
artifact root and keeps the method-specific layer and embedding keys explicit
in `result_inputs.tsv`.

## Inputs

`config.json` contains the dataset labels, cluster edges, BATC bin counts,
metric settings, and output names. `result_inputs.tsv` contains one row per
method/dataset pair. Its paths are relative to `--artifact-root`; it is the
runtime input map, not a history of result moves.

The runtime map is metric-specific and is intentionally separate from the
cell-render manifest in `manifests/cell_result_inputs.tsv`. For example, the
cell metrics use the DeepVelo metric-result family, while unified streamline
figures may use its `results_with_vu` family.

The evaluator copies the configured native velocity layer to an in-memory
`velocity` layer. It uses the configured source embedding. When the input
does not store a velocity embedding, such as the final DeepKINET results, the
embedding is recomputed from the native velocity layer after neighbors and
the velocity graph are prepared.

## Metrics

`metrics.py` calculates Consistency, CBDir, ICCoh, and BATC. Consistency uses
the mean cosine similarity between each cell and its 30 nearest neighbors.
CBDir and ICCoh use the UnitVelo evaluation implementation with the configured
cluster edges. BATC uses the dataset-specific bin count in `config.json`.

`transprob.py` calculates the companion TransProb/TransCosine score using the
same labels, neighbors, and velocity graph. Existing valid graphs are reused;
otherwise the graph is computed from the configured velocity layer.

`plots.py` consumes only the two compact result tables and renders the five
metric panels at the requested DPI. It does not reopen H5AD files.

## Run

Validate result paths and required input rows:

```bash
python evaluation/cell/metrics.py \
  --artifact-root /path/to/artifacts \
  --check-only
```

Write the four cell-level metrics:

```bash
python evaluation/cell/metrics.py \
  --artifact-root /path/to/artifacts \
  --output /path/to/metrics/all_metrics.csv
```

Write TransProb separately:

```bash
python evaluation/cell/transprob.py \
  --artifact-root /path/to/artifacts \
  --output /path/to/metrics/transprob.csv
```

Render the five panels:

```bash
python evaluation/cell/plots.py \
  --metrics /path/to/metrics/all_metrics.csv \
  --transprob /path/to/metrics/transprob.csv \
  --output /path/to/figures
```

The committed `results/` tables are the current final compact benchmark
values. Large H5AD files, checkpoints, logs, and historical audit manifests
remain external or under `provenance/`.

For the private curation contract check, use:

```bash
python manifests/validate_cell_contract.py
```
