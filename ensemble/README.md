# Signed-transition ensemble

Mean-4 and RowGeometricMedian-4 (GeoMedian-4) aggregate four selected Cell outputs:
TIvelo, veloVI without-vu, the DeepVelo root no-vu result family, and UniTVelo.
They are post-hoc comparators outside the 20-method configuration registry.

## Results and reproduction

| File | Use |
| --- | --- |
| [results/scores.csv](results/scores.csv) | Dataset-level scores for the two variants |
| [results/rankings_vs_20.csv](results/rankings_vs_20.csv) | Per-dataset ranks against the 20 Cell configurations |
| [results/metadata.json](results/metadata.json) | Protocol, input identities, frames, aggregate scores and run provenance |
| [results/summary.md](results/summary.md), [ranking_summary.md](results/ranking_summary.md) | Human-readable result summaries |

Run from the repository root with the inputs in [config.json](config.json).

Use the [ensemble dependency record](../docs/environments.md#evaluation).
`--artifact-root` is the velocity directory itself: each `results/...h5ad`
in the configuration is read below that root. The reference and four method
inputs can be located through the published member manifest and the
[local-root guide](../data/README.md#local-data-roots).

```bash
python ensemble/benchmark.py \
  --artifact-root /path/to/artifacts \
  --output-root outputs/ensemble \
  --dataset capofibro
```

Repeat `--dataset` to select a subset; omission runs all configured datasets.
`--save-streams` adds stream panels. `--h5ad-root` writes native H5ADs to an
external staging directory. Compact scores, metadata and summaries use the
output root. [plot_supfig10.py](plot_supfig10.py) renders S10 from the frozen
Cell/ensemble tables; [the figure map](../docs/figure_map.md) links final assembly.

Check dataset/variant coverage in `scores.csv` and the recorded input and
protocol identities in `metadata.json` under the selected output root.
For new method outputs, use a separate `--config` with their relative paths;
the default configuration selects the frozen inputs.

## Aggregation

The scVelo 0.3.4 post-processing path is:

```python
scv.tl.velocity_graph(adata, vkey="velocity", xkey="Ms")
scv.tl.transition_matrix(
    adata, scale=10, self_transitions=True, use_negative_cosines=True,
    weight_diffusion=0, n_neighbors=None,
)
```

The committed Cell reference graph supplies fixed 30-neighbor candidate edges
(`--reference-neighbor-root` selects another explicit copy).
Velocity-graph cosines use each method's native `Ms`/`velocity`.
Invalid or signed-L1-near-zero (≤1e-12) method rows become identity rows;
fallback counts are recorded. Off-diagonal operators align on a deterministic
union support.

- Mean-4 takes the arithmetic mean of each aligned edge.
- GeoMedian-4 takes the Euclidean geometric median of the four complete signed
  off-diagonal rows for each source cell.

Both are permutation-invariant and learn no method-identity weights.
Geometric-median convex coefficients are row-wise diagnostics.
Signed values and native directional mass are retained without post-hoc
normalization; remaining signed-L1 mass is put on the diagonal.
scVelo `velocity_embedding` receives the aggregate T explicitly.

## Coordinates and ranking

Schema-7 Forebrain outputs retain `X_pca`/`ensemble_velocity_pca` (1720 × 30)
for CBDir and `X_umap`/`ensemble_velocity_umap` (1720 × 2) for BATC.
The frozen Cell frames are checked by coordinates, dimensions and cell order;
no UMAP optimization occurs. Projection uses `autoscale=True`,
`all_comps=True`, `direct_pca_projection=False`.
Forebrain streams show PC1/PC2; the other datasets keep their UMAP spaces.
The induced gene-space velocity is a projection adapter.

For each dataset and metric, rank is
`1 + count(Cell method score > ensemble score)`, with competition ties and
`tie_count`. Each variant is compared separately with the 20 Cell methods
(21 participants); variants do not enter each other's candidate set.
Ranking is dataset-specific, not a cross-dataset average.
The config alias `dentate_lamanno` maps to Cell `DentateGyrus_lamanno`.

Historical schema-6 Forebrain fields used PCA under UMAP names.
Their original bytes remain retained; current schema-7 artifacts distinguish
both frames. Numerical and migration evidence is kept in the
[registered historical record](../docs/fixed_neighbor_reproduction_20260912.md)
and current metadata, outside the instructions needed for a new run.
