# Cell-level velocity evaluation

Evaluate cell-level velocity predictions from 20 configurations on 14 real
datasets using five metrics.
The committed [primary table](results/cell_metrics.csv) has 280 rows.
[Simulation Cell](simulation/README.md) and the
[gene-output-matched comparison](#gene-output-matched-cell-comparison) are
separate analyses.

## Inspect and plot results

| File | Contents |
| --- | --- |
| [results/cell_metrics.csv](results/cell_metrics.csv) | Consistency, ICCoh, CBDir, TransCosine and BATC |
| [results/cell_metrics_supplementary.csv](results/cell_metrics_supplementary.csv) | Additional metrics, CBDir-HD (`highdim_cbdir`) and evaluability counts |
| [results/cell_metrics_audit.csv](results/cell_metrics_audit.csv) | Source identities, layer keys, invalid-value counts and status |
| [results/cell_metrics_edge_details.csv](results/cell_metrics_edge_details.csv) | Transition-level results |
| [results/tivelo_parity.csv](results/tivelo_parity.csv) | Comparison with the official TIvelo scorer |
| [streams/](streams/) | Retained stream images and source mapping |

Run from the repository root:

```bash
python evaluation/cell/validate_results.py
python evaluation/cell/plots.py --output outputs/cell-plots
```

The validator checks the committed bundle. The plotter consumes its compact
tables; `--metrics` and `--supplementary-metrics` select other result tables.

## Recalculate metrics

Materialize the [saved method outputs](../../data/RESULT_DOWNLOADS.md) and
reference moments. `--artifact-root` points to the parent of `velocity/`;
each `velocity:<path>` in [inputs.tsv](inputs.tsv) resolves to
`<artifact-root>/velocity/<path>`. The reference-moments directory contains
`<evaluation_dataset>.npz`; neighbors and shared embedding assets default
to their committed directories.

Select two local worker interpreters using the
[evaluation dependency records](../../docs/environments.md#evaluation):

| Option | Work performed |
| --- | --- |
| `--indicators-python` | Runs `prepare_cell_metrics.py`: aligns inputs, installs the fixed reference graph, prepares expression/velocity arrays, computes Consistency/ICCoh/CBDir/BATC and the local TIvelo parity values |
| `--tivelo-python` | Runs `score_tivelo_metrics.py` on the prepared H5AD using TIvelo 0.1.4's official TransCosine and CBDir-HD functions |

Supply executable paths from your installation. Environment names are arbitrary:

```bash
CELL_METRICS_PYTHON=/path/to/cell-evaluation/bin/python
TIVELO_METRICS_PYTHON=/path/to/tivelo-metrics/bin/python
```

The controller can run with Python 3.9+ and invokes these workers in sequence:

```bash
python evaluation/cell/run_cell_metrics.py \
  --artifact-root /path/to/artifacts \
  --reference-moments-root /path/to/reference-moments \
  --output-root outputs/cell-evaluation \
  --indicators-python "$CELL_METRICS_PYTHON" \
  --tivelo-python "$TIVELO_METRICS_PYTHON" \
  --prepare-script evaluation/cell/prepare_cell_metrics.py \
  --official-script evaluation/cell/score_tivelo_metrics.py
```

[inputs.tsv](inputs.tsv) and [config.json](config.json) define the selected
files, fields, methods, datasets and transitions.
For MultiVelo, [native_outputs.tsv](native_outputs.tsv) binds each Cell output
to its native fit. Cell scoring selects `velo_s_norm`, `velo_s_norm_genes`,
`velo_s_norm_params` and `velo_s_norm_umap` from the official MultiVelo graph
and projection workflow. The global normalization follows MultiVelo 0.1.0.
The native fit files used for Gene scoring remain separately identified.
Use `--inputs` / `--config` for other scopes and a separate output root for
partial runs. The canonical Cell directory requires all 280 pairs.
`--method` uses operational IDs and `--dataset` uses `evaluation_dataset`
values from `inputs.tsv`; both options can be repeated. When a dataset's
shared frame is read from scVelo deterministic, that reference H5AD is also
needed for a single-method run. Keep its row in a copied input map.

For a new method result, copy the input map, update its `result_path`,
`sha256`, size and consumed layer fields, and pass the copy with `--inputs`.
The logical path must resolve below the chosen artifact root. This keeps the
new result distinct from the frozen benchmark input.

`--reference-neighbors-root` defaults to the committed
[reference_neighbors/](reference_neighbors/).
Temporary pair files and logs use `--work-root`, or an ignored `.work`
directory under the output root.

The runner writes primary/supplementary scores, source audits, edge details and
TIvelo parity. Extended edge details retain native CBDir edge means/counts and
BATC pre-max candidate means/winner counts. Complete Forebrain runs also write
its PCA edge table and manifest. Resume requires matching inputs, code,
configuration, metric spaces and extended details.

Results are written to `<output-root>/results/`, including
`cell_metrics.csv`, `cell_metrics_supplementary.csv`, `cell_metrics_audit.csv`,
`cell_metrics_edge_details.csv` and `tivelo_parity.csv`. Check the selected
dataset/method rows, source identities and parity outcomes there; per-pair
logs are below the work root. `validate_results.py` checks the committed
canonical bundle. To plot another result table, use the plotter's
`--metrics` and `--supplementary-metrics` input options.

## Metric coordinates and numerical policy

The stored reference graphs use 30 neighbors, 30 PCs, Euclidean distance,
`X_pca` and random state 0. Installing a graph transfers its neighbor indices,
distances and connectivities. For the primary Cell benchmark, velocity-graph
cosines use reference `Ms` aligned to each method's cells and genes, together
with its prepared velocity. veloAE retains its native latent state and velocity.

Forebrain CBDir uses the full shared PCA30 velocity projection and streamlines
show its first two components. BATC uses the frozen shared UMAP2.
[reference_embeddings/](reference_embeddings/) and
`reference_neighbors/forebrain.npz` retain those frames with cell-order and
coordinate checks. Exact-zero and empty CBDir edges are excluded.

BATC follows GRAVITY's curve fitting, tangent-cosine and per-cell branch
maximization. NaN branch scores become zero before maximization.

Stored CSVs retain full precision. Individual/summary annotations and formal
tables use four decimals, dense heatmaps two; ranking uses unrounded scores.

## Gene-output-matched cell comparison

[gene_matched/inputs.tsv](gene_matched/inputs.tsv) selects Gene-analysis outputs
for DeepVelo, UniTVelo, veloVI and MultiVelo on seven real datasets: 28 scored rows, plus
seven reference-only scVelo locators for the canonical UMAP. It uses the same
alignment, fixed graph, reference moments, metrics and parity checks as Cell.

- DeepVelo changes `pred_unspliced=false` (Cell) to `true` (Gene/Kinetics),
  changing the training head and spliced output.
- UniTVelo adds four Gene layers; shared arrays and all six Cell metrics match.
- MultiVelo uses the official normalization of native `velo_s`, selecting
  `velo_s_norm` with its corresponding gene mask, parameters and projection.
  Primary and matched Cell evaluations use this same output on the seven
  shared datasets; Gene scoring uses native `velo_s` and `velo_u`.
- veloVI Cell and Gene are distinct training runs. Both register native
  `Ms`/`Mu`. Cell uses package-default training (up to 500 epochs, batch 256)
  and divides spliced velocity by `20 / max(latent_time)` per gene.
  Gene uses up to 1,000 epochs, batch 64 GPU/32 CPU, unscaled spliced velocity
  and the trained model's unspliced extraction.

These comparisons evaluate the complete output workflows described above.
The 28 supplementary rows are stored separately from the primary 280-row table.
Together with the other five methods' primary Cell rows, they reproduce S9's
315 scores (nine methods × seven datasets × five displayed metrics).
Render its five Cell panels directly from these current tables:

```bash
python evaluation/cell/plots.py --gene-matched --output outputs/s9-cell-panels
```


```bash
python evaluation/cell/run_cell_metrics.py \
  --artifact-root /path/to/artifacts \
  --inputs evaluation/cell/gene_matched/inputs.tsv \
  --reference-moments-root /path/to/reference-moments \
  --output-root outputs/cell-gene-matched \
  --indicators-python "$CELL_METRICS_PYTHON" \
  --tivelo-python "$TIVELO_METRICS_PYTHON" \
  --prepare-script evaluation/cell/prepare_cell_metrics.py \
  --official-script evaluation/cell/score_tivelo_metrics.py \
  --method deepvelo --method unitvelo --method veloVI --method multivelo
```

## Website previews

[streams/preview_manifest.json](streams/preview_manifest.json) binds the 280
previews to exact input identities and their stored Cell scores.
Rebuild in a separate checkout using the
[Cell preparation dependencies](../../docs/environments.md#evaluation), including scVelo 0.2.5:

```bash
python evaluation/cell/render_stream_previews.py \
  --artifact-root /path/to/artifacts \
  --reference-moments-root /path/to/reference-moments \
  --cache-root site/tmp/stream-frame-cache \
  --manifest-root site/tmp/stream-manifests
```

Use the frozen [palette](stream_preview_config.json) and assemble the complete
dataset fragments before [website export](../../site/README.md).
Figure 2 and Supplementary Figure 12 stream rows can be rendered from
[Source Data](../../data/SOURCE_DATA.md) display arrays, with the selected
MultiVelo native projection supplied by `inputs.tsv`:

```bash
python evaluation/cell/plot_cell_stream_panels.py \
  --source-data-root /path/to/source-data \
  --artifact-root /path/to/velocity-artifacts \
  --output-root outputs/cell-stream-panels
```

`--dataset` selects a registered stream panel; omission renders all five.
The two real-data reference panels retain their annotated prior geometry in
[figures/references/](figures/references/).
Forebrain uses the reference-30 projection; other previews reuse the selected
native velocities and reference coordinates.
