# Normalized five-member ensemble

The selected **GeoMedian-5** combines TIvelo, veloVI, DeepVelo, UniTVelo and
GraphVelo_scv on the 14 benchmark datasets. The first four members supply
signed transition operators `T`, and GraphVelo_scv supplies its native
coefficient matrix `Phi`. The replay uses veloVI without-vu and the DeepVelo
Cell no-vu result family. See [config.json](config.json), the
[method registry](../methods/configurations.tsv) and the
[dataset registry](../data/datasets.tsv) for input identities.

## Selected results

| Artifact | Purpose |
| --- | --- |
| [scores.csv](results/selected/scores.csv) and [summary.csv](results/selected/summary.csv) | Selected per-dataset CBDir/BATC and equal-dataset macro means and medians |
| [Table S5 A](results/selected/Table_S5_A.csv) | Full ensemble, TIvelo + veloVI, five leave-one-out configurations and Mean-5, with paired deltas and Pareto membership |
| [Table S5 B](results/selected/Table_S5_B.csv) | Selected CBDir/BATC for each dataset and the macro-average |
| [dataset_scores.csv](results/selected/dataset_scores.csv) | All 112 dataset-configuration scores and evaluation counts |
| [projections](results/selected/projections/) | Selected projected fields, cell IDs, labels, frozen neighbors and coordinates |
| [protocol.json](results/selected/protocol.json) and [metadata.json](results/selected/metadata.json) | Settings, versions, input and code identities, and output hashes |

The projection arrays retain `ensemble_velocity_umap` and, for Forebrain,
`ensemble_velocity_pca`. They represent velocity fields in the evaluation
coordinate frames. Selection and evaluation use the same 14 datasets, so the
comparisons are exploratory.

## Scientific computation

1. Verify the frozen input manifest, geometry and five signed sparse
   operators per dataset. The four saved `T` inputs were constructed from
   native expression and velocity state with the fixed neighbor graph.
2. Discard diagonal coefficients and divide each active row by its
   off-diagonal absolute weight sum, preserving coefficient signs. Rows
   with mass at or below `1e-12` contribute zero.
3. Align the signed rows for each source cell on the union of their
   nonzero supports. Take the equal-weight Euclidean geometric median,
   or the arithmetic mean for Mean-5.
4. Apply GraphVelo's weighted raw displacement
   `V_i = sum_j W_ij * (Y_j - Y_i)` without centering, diagonal completion
   or post-aggregation normalization.
5. Evaluate CBDir and BATC using the fixed neighbors and annotations.
   Forebrain CBDir uses PCA30 and BATC uses UMAP2. Both metrics use UMAP2
   for the other datasets.

## Reproduce the selected bundle

Use [the ensemble evaluation environment](../docs/environments.md#evaluation).
Run from the repository root and choose a new output directory.

```bash
PYTHONNOUSERSITE=1 python -m ensemble --jobs 3 --output-root outputs/geomedian5
```

The replay reads `ensemble/frozen_inputs/` by default. Use `--input-root`
to point to another copy of the frozen bundle, which must include
`manifest.json` and one geometry archive and five sparse operator archives
per dataset. Download [the frozen five-member inputs](https://huggingface.co/datasets/Zelin259/rna-velocity-benchmark/resolve/v1.1.0/v1.1.0/ensemble_outputs/geomedian5_frozen_inputs.tar.gz?download=true), extract the archive into a new directory, and use that directory as `--input-root`. The [result-package manifest](../data/result_package_manifest.tsv) records its SHA-256. The command consumes
local files and performs no downloads or method training. It writes the selected results, Table S5 and the Figure
S10 source tables. The replay produces tabular data and projected arrays.

## Historical four-member workflows

The [archived GeoMedian-4 bundle](results/legacy_geomedian4_graphvelo/)
used TIvelo, veloVI, DeepVelo and UniTVelo with their native off-diagonal
row mass and GraphVelo's weighted-displacement projection. Its retained
[stream atlas](results/legacy_geomedian4_graphvelo/streams.pdf) includes
both Forebrain views. Table S5's matched `w/o GraphVelo_scv` comparison
uses the same four member identities with the selected per-member signed-L1
normalization. The [legacy scVelo-projection results](results/scores.csv)
are a separate four-member baseline.

For workflows
that depend on rebuilding those legacy H5AD operators from saved method
outputs, use [benchmark.py](benchmark.py) with `--artifact-root`
(or opt into the historical builder with `--legacy-four`):

```bash
python ensemble/benchmark.py \
  --artifact-root /path/to/velocity \
  --output-root outputs/four-method-operators \
  --h5ad-root /path/to/new-four-member-h5ads \
  --n-jobs 8
```

The builder writes the historical scVelo projection and both aggregation
controls. Newly generated H5ADs require their own matching release manifest
before projection replay.

## Dataset resampling

[pareto_bootstrap.py](pareto_bootstrap.py) uses the formal 20-configuration Cell score table. In each of 10,000 draws (seed 20261109), all configurations and both directional metrics share a sample of 14 datasets. [S10a](results/selected/s10/panel_a_data.csv) records the Pareto-front inclusion frequencies. The selected member and aggregation controls are in [S5A](results/selected/Table_S5_A.csv).

## Supplementary Figure 10

The [replay tables](results/selected/s10/) contain the data for all three
panels. Panel A uses 10,000 paired dataset bootstrap draws with seed
`20261109` to estimate Pareto-front frequencies for the 20 single-method
configurations. Panels B/C compare their CBDir/BATC with GeoMedian-5 on
14 datasets. Each score panel orders methods by its own macro mean and
marks the five ensemble members with stars.

The [editable assembly](../figures/assembly/suppfig10/suppfig10.pptx) and
its [panel manifest](../figures/assembly/suppfig10/panel_manifest.json) record
the approved layout and plotted data. The [current author PDF](../figures/suppFigTen.pdf) was exported from that assembly. The standalone metric renderer accepts the historical
four-member score table. To regenerate those historical panels:

```bash
python ensemble/plot_supfig10.py \
  --ensemble-csv ensemble/results/legacy_geomedian4_graphvelo/scores.csv \
  --output-dir outputs/suppfig10-geomedian4
```
