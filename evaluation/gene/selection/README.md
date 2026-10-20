# Real-data gene selection

## Final selection and evidence

[selected_genes.tsv](selected_genes.tsv) identifies the 59 frozen dataset–gene
pairs evaluated by [the Gene entry point](../run_reference30_gene_cbdir.py).
Each pair has candidate-table support. The final genes were selected manually
after slope screening and phase-portrait inspection. Row ranks in the TSV
locate the supporting candidate rows. Synthetic genes
follow the simulator's known modules.

| Material | Role |
| --- | --- |
| [selected_genes.tsv](selected_genes.tsv) | Final evaluated pairs and candidate evidence |
| [provenance.json](provenance.json) | Original source locators, script snapshots, checksums and limitations |

## Screening rules

### MURK

[Barile et al. (2021)](https://doi.org/10.1186/s13059-021-02414-y) provides the
[upstream notebook](https://github.com/mebarile/Gata1_Erythroid_kinetics),
pinned to commit `1e0de6338a85a559478627927c24bdc72d2d8919`; the unmodified
notebook and GPL-3.0 license remain in [upstream/](upstream/).

The notebook combines Scanpy variable genes with scVelo-filtered genes,
computes spliced/unspliced moments, fits population-specific slopes, applies
expression and slope-change criteria, and ranks scaled slope differences.
It does not directly estimate transcription-rate alpha. Its interval uses
`t.ppf(0.95, n-1)`, unlike the branching two-sided interval.

[reproduce_murk.py](reproduce_murk.py) executes the unchanged numerical cells
with the input read redirected.

### Branching

[reproduce_branching.py](reproduce_branching.py) adapts the slope-comparison
idea, combining the two variable-gene sets and scaling moments by the maximum
across both populations. It fits OLS slopes with an intercept and two-sided
95% Student-t intervals (`t.ppf(0.975, n-2)`), then ranks non-overlapping
interval gaps without a positive-gap cutoff.

The `pro` filter requires maximum absolute expression-versus-cluster-rank
Spearman correlation strictly greater than 0.1 as an expression-trend magnitude
threshold. Retina retains the historical
equal Neuroblast/RGC ranks, producing zero correlation in that branch.

| DentateGyrus task | Populations and ranks |
| --- | --- |
| Neuron/Glia | Neuron: `RadialGlia, nIPC, Nbl1, Nbl2, ImmGranule1, ImmGranule2, Granule, CA, CA1-Sub, CA2-3-4`; Glia: `RadialGlia, GlialProg, ImmAstro, OPC` |
| Granule/Pyramidal | `Nbl2, ImmGranule1, ImmGranule2, Granule` versus `Nbl2, CA, CA1-Sub, CA2-3-4`; ranks 1/2/3/4 versus 1/2/3/3 |

Selection-stage `moments(n_neighbors=30)` is distinct from the fixed-reference-30
evaluation graph.

## Reconstruct candidate tables

Run from the repository root with the input H5ADs and compatible
Scanpy/scVelo environments. Manifests record package versions, input hashes,
hash seed and preprocessing settings; the retained two-task run's thread
settings are indexed in the summary. Historical `list(set(...))` ordering
depends on Python's hash seed. Each screening entry requires a new output directory:

```bash
python evaluation/gene/selection/reproduce_murk.py \
  --input /path/to/Gastrulation.h5ad \
  --input-id artifact://velocity/data/Gastrulation/Gastrulation.h5ad \
  --output /path/to/new-murk-check

PYTHONHASHSEED=0 python evaluation/gene/selection/reproduce_branching.py \
  --dataset DentateGyrus_lamanno --include-neuron-glia \
  --input /path/to/DentateGyrus_lamanno_for_scvelo.h5ad \
  --input-id artifact://velocity/data/DentateGyrus_lamanno/DentateGyrus_lamanno_for_scvelo.h5ad \
  --output /path/to/new-dentategyrus-check

PYTHONHASHSEED=0 python evaluation/gene/selection/reproduce_branching.py \
  --dataset retina --input /path/to/retina.h5ad \
  --input-id artifact://velocity/indicators/gene/Mouse_retina_analyze/data/retina.h5ad \
  --output /path/to/new-retina-check
```

`--include-neuron-glia` adds the recovered December Neuron/Glia task to
Granule/Pyramidal; its filtered output uses the January rule. Without the
option, only Granule/Pyramidal outputs are emitted. Retina emits unfiltered
and filtered PR/RGC tables.

## Historical agreement

| Material | Role |
| --- | --- |
| [historical/](historical/) | Original candidate tables and source snapshots |
| [results/summary.json](results/summary.json) | Reconstruction attempts and full-precision comparisons |
| [historical/MURK_selected.md](historical/MURK_selected.md) | Broader intermediate lists; the TSV and evaluation entry define the final selection |

Per-gene selection decision logs are unavailable. Recovered MURK population
comparisons apply to Gastrulation; configurations for the other four MURK
datasets remain missing. The legacy Gastrulation input is absent, so the
retained raw-input run is a reconstruction check rather than exact
historical-input reproduction. The commented project excerpt is configuration
evidence, not a run log.

The [preserved snapshots](historical/scripts/) separate the December 2025
unfiltered tasks, January 7 filtered Granule/Pyramidal ranks, and January 18
ten-gene list. The two slope tables record separate population comparisons;
the later list records the manual selection.

| Check | Agreement with historical tables |
| --- | --- |
| Retina branching | All 2,264 candidates and saved numerical fields match exactly; primary-score ranks match after accounting for ties. CSV ordering differs only within equal-score groups. |
| DentateGyrus branching | Both slope tasks recover all 3,841 candidates; the filtered table recovers all 3,300 candidates. Each top-50 set overlaps 50/50. Values and non-tied ranks are not exact. All ten final genes occur in all three tables. |
| Gastrulation MURK | The run returns 89 candidates containing all 45 historical candidates. The historical list is a manually selected subset, not a top-45 rule. Shared numerical values are not exact. |
| Other MURK datasets | Historical tables are retained; dataset-specific population configurations were not recovered, so screening was not rerun. |

Top-50 overlap reports set membership. DentateGyrus maximum interval-gap
discrepancies are approximately 0.0245 (Neuron/Glia), 0.0321
(Granule/Pyramidal), and 0.0322 (filtered). All 59 selected pairs have candidate
support; the table above records the numerical agreement of each reconstruction.

[compare_rankings.py](compare_rankings.py) separates membership, row order,
tie-aware ranks, top-50 overlap and full-precision differences. Use a new
output file for comparison:

```bash
python evaluation/gene/selection/compare_rankings.py \
  --historical /path/to/historical.csv --rerun /path/to/reconstructed.csv \
  --output /path/to/new-comparison.json
```
