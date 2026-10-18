# Cell-level evaluation

Evaluate five trajectory metrics for 20 configurations on 14 real datasets.
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
reference moments, and provide the indicators and TIvelo Python environments:

```bash
python evaluation/cell/run_cell_metrics.py \
  --artifact-root /path/to/artifacts \
  --reference-moments-root /path/to/reference-moments \
  --output-root outputs/cell-evaluation \
  --indicators-python /path/to/indicators-python \
  --tivelo-python /path/to/tivelo-python \
  --prepare-script evaluation/cell/prepare_cell_metrics.py \
  --official-script evaluation/cell/score_tivelo_metrics.py
```

[inputs.tsv](inputs.tsv) and [config.json](config.json) define the selected
files, fields, methods, datasets and transitions.
Use `--inputs` / `--config` for other scopes and a separate output root for
partial runs. The canonical Cell directory requires all 280 pairs.
`--reference-neighbors-root` defaults to the committed
[reference_neighbors/](reference_neighbors/).
Temporary pair files and logs use `--work-root`, or an ignored `.work`
directory under the output root.

The runner writes primary/supplementary scores, source audits, edge details and
TIvelo parity. Extended edge details retain native CBDir edge means/counts and
BATC pre-max candidate means/winner counts. Complete Forebrain runs also write
its PCA edge table and manifest. Resume requires matching inputs, code,
configuration, metric spaces and extended details; old incompatible caches
are recomputed. Incomplete Forebrain staging runs remove stale complete-run
audit claims.

## Metric coordinates and numerical policy

The stored reference graphs use 30 neighbors, 30 PCs, Euclidean distance,
`X_pca` and random state 0. Installing a graph transfers its neighbor indices,
distances and connectivities. Method-native expression/velocity spaces supply
the velocity-graph cosines.

Forebrain CBDir uses the full shared PCA30 velocity projection and streamlines
show its first two components. BATC uses the frozen shared UMAP2.
[reference_embeddings/](reference_embeddings/) and
`reference_neighbors/forebrain.npz` retain those frames with cell-order and
coordinate checks. Native H5ADs are unchanged. The 120 Forebrain transition
rows retain exact-zero/empty-edge exclusion; DeepVelo Cell uses its no-vu output.

BATC follows GRAVITY's curve fitting, tangent-cosine and per-cell branch
maximization, with one benchmark adaptation: NaN branch scores become zero
before maximization. A zero fallback can exceed a finite negative score.
Upstream `np.nanargmax` is undefined for all-NaN sources in
`bonemarrow/TSvelo` and `rpe1_gen/TFvelo`; the other 278 pairs matched the
benchmark exactly under the original-policy check, as did all 28 ensemble
scores. See the [historical numerical record](../../docs/fixed_neighbor_reproduction_20260912.md).

Stored CSVs retain full precision. Individual/summary annotations and formal
tables use four decimals, dense heatmaps two; ranking uses unrounded scores.

## Gene-output-matched cell comparison

[gene_matched/inputs.tsv](gene_matched/inputs.tsv) selects Gene-analysis outputs
for DeepVelo, UniTVelo and veloVI on seven real datasets: 21 scored rows, plus
seven reference-only scVelo locators for the canonical UMAP. It uses the same
alignment, fixed graph, reference moments, metrics and parity checks as Cell.

- DeepVelo changes `pred_unspliced=false` (Cell) to `true` (Gene/Kinetics),
  changing the training head and spliced output.
- UniTVelo adds four Gene layers; shared arrays and all six Cell metrics match.
- veloVI Cell and Gene are distinct training runs. Both register native
  `Ms`/`Mu`. Cell uses package-default training (up to 500 epochs, batch 256)
  and divides spliced velocity by `20 / max(latent_time)` per gene.
  Gene uses up to 1,000 epochs, batch 64 GPU/32 CPU, unscaled spliced velocity
  and the trained model's unspliced extraction. Extraction is not a training
  switch. Dentategyrus_lamanno outputs share 1,290 of their 1,291 genes.

These comparisons measure output-workflow sensitivity; they do not isolate
appending an unspliced layer or scaling alone.
The supplementary results leave the canonical 280 rows unchanged.

```bash
python evaluation/cell/run_cell_metrics.py \
  --artifact-root /path/to/artifacts \
  --inputs evaluation/cell/gene_matched/inputs.tsv \
  --reference-moments-root /path/to/reference-moments \
  --output-root outputs/cell-gene-matched \
  --indicators-python /path/to/indicators-python \
  --tivelo-python /path/to/tivelo-python \
  --prepare-script evaluation/cell/prepare_cell_metrics.py \
  --official-script evaluation/cell/score_tivelo_metrics.py \
  --method deepvelo --method unitvelo --method veloVI
```

## Website previews

[streams/preview_manifest.json](streams/preview_manifest.json) binds the 280
previews to exact input identities. Previews add displays, not scores.
Rebuild with the scVelo 0.2.5 indicators environment in a separate checkout:

```bash
python evaluation/cell/render_stream_previews.py \
  --artifact-root /path/to/artifacts \
  --reference-moments-root /path/to/reference-moments \
  --cache-root site/tmp/stream-frame-cache \
  --manifest-root site/tmp/stream-manifests
```

Use the frozen [palette](stream_preview_config.json) and assemble the complete
dataset fragments before [website export](../../site/README.md).
Forebrain uses the reference-30 projection; other previews reuse the selected
native velocities and reference coordinates.
