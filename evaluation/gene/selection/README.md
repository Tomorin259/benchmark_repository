# Real-data gene selection

This module preserves computational candidate screening and the subsequent
manual gene selection used for the real-data gene-level benchmark. Screening
outputs do not replace the manuscript gene lists or reference-30 velocity metrics.

## Sources and final selection

- `provenance.json` records original source locators, dated script snapshots,
  checksums, and computational limitations.
- `historical/` preserves original candidate tables and source snapshots.
  Snapshots retain their historical paths and are evidence, not public launchers.
- `selected_genes.tsv` records the 59 dataset-gene pairs evaluated by
  `../run_reference30_gene_cbdir.py`. Every pair occurs in a recovered candidate
  table. Its first-evidence row rank is not a combined ranking or top-k rule.
- `results/` contains the retained checks, input/environment manifests,
  and full-precision comparisons against the historical tables.
  `results/summary.json` indexes these checks; each comparison JSON is
  the numerical record for its associated historical/rerun table pair.

Results use stable gene-setting and dataset names:

```text
results/
├── branching/
│   ├── DentateGyrus_lamanno/
│   └── retina/
├── murk/
│   └── gastrulation/
└── summary.json
```

Package versions, hash seed, input checksums, and run settings belong in the
manifests, not in directory names. These remain reconstruction checks rather
than replacement benchmark results.

The authors selected the final genes manually after computational screening
and phase-portrait inspection. Per-gene decision logs were not retained.
No automated intersection, union, or combined-score rule is claimed for this
final step. The historical `MURK_selected.md` contains broader intermediate
lists; the current evaluation entry and `selected_genes.tsv` define the final
selection. Synthetic genes instead follow the simulator's known modules.

## MURK screening

The upstream source is Barile et al. (2021),
https://doi.org/10.1186/s13059-021-02414-y, with the notebook pinned to commit
`1e0de6338a85a559478627927c24bdc72d2d8919` of
https://github.com/mebarile/Gata1_Erythroid_kinetics.
The unmodified notebook and GPL-3.0 license are retained in `upstream/`.

The notebook combines Scanpy variable genes and scVelo-filtered genes, computes
spliced/unspliced moments, fits population-specific slopes, applies expression
and slope-change criteria, and ranks candidates by scaled slope differences.
It does not directly estimate transcription-rate alpha for selection.
Its uncertainty calculation uses `t.ppf(0.95, n-1)`, distinct from the
branching script's two-sided 95% interval.

`reproduce_murk.py` executes the unchanged numerical notebook cells, redirecting
only the input read. Its population comparisons apply to Gastrulation.
The comparisons used for the other four MURK datasets have not been recovered;
their historical candidate tables are preserved without applying guessed groups.
The legacy Gastrulation input is absent, so the retained run uses the currently
available raw input and is a reconstruction check, not exact historical-input
reproduction. The recovered commented project excerpt is configuration evidence,
not a run log.

```bash
python evaluation/gene/selection/reproduce_murk.py \
  --input /path/to/Gastrulation.h5ad \
  --input-id artifact://velocity/data/Gastrulation/Gastrulation.h5ad \
  --output /path/to/new-murk-check
```

## Branching screening

The project scripts adapt the phase-portrait slope-comparison idea rather than
the complete MURK rule. They combine Scanpy variable genes and scVelo-filtered
genes, scale moments by the maximum across both comparison populations, fit
OLS slopes with an intercept and two-sided 95% Student-t intervals
(`t.ppf(0.975, n-2)`), and rank by non-overlapping interval gaps.
There is no positive-gap cutoff.

The `pro` variants retain genes whose maximum absolute expression-versus-cluster-
rank Spearman correlation is strictly greater than 0.1. This is a magnitude
threshold, not a p-value test; opposite-signed trends are not required.
The Retina snapshot assigns the same rank to Neuroblast and RGC, yielding zero
correlation in that branch. This historical configuration is retained.

Three dated DentateGyrus snapshots establish the recovered sequence:

- `branching_DentateGyrus_lamanno_20251216.py` contains the unfiltered
  Neuron/Glia and Granule/Pyramidal slope tasks and the `for_scvelo` input path.
- `branching_DentateGyrus_lamanno_pro_20260107.py` supplies the matching
  Granule/Pyramidal trend ranks and filtered-output filename.
- `branching_DentateGyrus_lamanno_analyze_20260118.py` already contains the
  current ten-gene list. It does not implement an automated final-selection rule.

For Neuron/Glia, Neuron contains `RadialGlia, nIPC, Nbl1, Nbl2, ImmGranule1,
ImmGranule2, Granule, CA, CA1-Sub, CA2-3-4`; Glia contains
`RadialGlia, GlialProg, ImmAstro, OPC`. Granule/Pyramidal compares
`Nbl2, ImmGranule1, ImmGranule2, Granule` against
`Nbl2, CA, CA1-Sub, CA2-3-4`, with ranks 1/2/3/4 and 1/2/3/3.
The two historical slope tables are separate comparisons, not a combined ranking.

```bash
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

The DentateGyrus entry emits both December slope-task outputs when
`--include-neuron-glia` is enabled, plus the Granule/Pyramidal filtered output
using the January rule. The latter is not attributed to the December script.
Without this option it emits only the Granule/Pyramidal outputs.
Retina emits unfiltered and filtered PR/RGC outputs.

Selection-stage `moments(n_neighbors=30)` is not the fixed-reference-30
evaluation graph. The historical `list(set(...))` gene union has uncontrolled
ordering unless Python's hash seed is fixed. Run manifests record the actual
hash seed, packages, input checksum, groups, and parameters; launch thread
settings for the retained two-task run are recorded in the summary.

## Retained results and verification

| Check | Agreement with historical tables |
| --- | --- |
| Retina branching | All 2,264 candidates and saved numerical fields match exactly; primary-score ranks match after accounting for ties. CSV ordering differs only within equal-score groups. |
| DentateGyrus branching | Both slope tasks recover all 3,841 candidates; the filtered table recovers all 3,300 candidates. Each top-50 set overlaps 50/50. Values and non-tied ranks are not exact. All ten final genes occur in all three tables. |
| Gastrulation MURK | The run returns 89 candidates containing all 45 historical candidates. The historical list is a manually selected subset, not a top-45 rule. Shared numerical values are not exact. |
| Other MURK datasets | Historical tables are retained; dataset-specific population configurations were not recovered, so screening was not rerun. |

Top-50 overlap denotes set membership, not exact ordering.
DentateGyrus maximum interval-gap discrepancies are approximately 0.0245 for
Neuron/Glia, 0.0321 for Granule/Pyramidal, and 0.0322 after trend filtering.
These checks establish candidate-table support for all 59 frozen pairs without
claiming exact historical numerical reproduction where discrepancies remain.
No rerun output is used to reselect genes or replace benchmark metrics.

All screening entry points require a new output directory.
`compare_rankings.py` reports candidate membership, row order, tie-aware ranks,
top-50 overlap, and full-precision numerical differences separately.
Its diagnostic tolerances do not convert discrepancies into exact agreement.

```bash
python evaluation/gene/selection/compare_rankings.py \
  --historical /path/to/historical.csv --rerun /path/to/reconstructed.csv \
  --output /path/to/new-comparison.json

python -m unittest discover -s tests -p test_gene_selection_provenance.py -v
```
