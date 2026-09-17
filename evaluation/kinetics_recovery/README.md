# Kinetics recovery and Figure 5

This module evaluates saved method-result H5AD files against the published
scEU-seq references for RPE1 and intestinal organoid. It does not prepare
datasets or train methods.

## Input policy and implementation status

The selected result source is the central `<velocity-root>/results/` directory.
Dynamo(label data) is the explicit exception: its labelled RPE1 result must
not be replaced by the label-free Dynamo(RNA_seq) H5AD. Freeze the selected
H5AD files and record their resolved paths and SHA256 before evaluation.
Evaluation and plotting must not overwrite or retrain these inputs.

The current evaluation is in `results/central_results/`; the canonical panels
and exact plotted rows are in `figures/fig5/`. The run manifest records selected
H5AD/reference hashes, software versions, renderer code hashes, and output
checksums. No methods were trained and `metrics.py` is unchanged.

Quantitative panels read one common evaluator output rather than independent
static legacy tables. Streams use saved embeddings and velocity projections
without recomputing neighbors, moments, graphs, or embeddings. Some saved
method files lack independent projection provenance; this is disclosed in the
stream manifests. Dynamo A is the retained labelled-data manual exception.

`../../figures/figFive.pdf`, assembly copies, assembly hashes, and the public
artifact registry are synchronized with the current panels. The composite is
a reviewed PowerPoint export with JPEG-encoded panels, not byte-identical PNG
embedding. `../../figures/fig5_composite_manifest.json` records its checksum,
panel mapping, and image-content validation. Old panels, the old composite,
and comparison materials are kept in an external saved archive, not alongside
current figures. `results/final_seed123/` is historical provenance retained for
the separately frozen Figure 6 input tables; it is not the current Figure 5
evaluation source. Figure 6 was not regenerated in this update.

## Selected H5AD files and panel coverage

Paths below are relative to the user-supplied velocity root. The external
kinetics archive can expose these through `method-results/kinetics/` aliases;
the resolved H5AD and hash, rather than the alias name, define the input.

| Method | System / purpose | Selected H5AD | Values consumed | Panels |
| --- | --- | --- | --- | --- |
| cellDancer | RPE1 | `results/celldancer/celldancer_rpe1.h5ad` | `celldancer_alpha`, `celldancer_gamma`; stream `velocity` | A, B, C, D |
| cellDancer | organoid | `results/celldancer/celldancer_organoid.h5ad` | Same layer keys | E, F |
| DeepVelo | RPE1 | `results/deepvelo/results_with_vu/deepvelo_rpe1.h5ad` | `cell_specific_alpha`, `cell_specific_gamma`; stream `velocity` | A, B, C, D |
| DeepVelo | organoid | `results/deepvelo/results_with_vu/deepvelo_organoid.h5ad` | Same layer keys | E, F |
| DeepKINET | RPE1 | `results/DeepKINET/DeepKINET_rpe1.h5ad` | `degradation_rate` as gamma; stream `DeepKINET_velocity` | A, B gamma, C gamma |
| DeepKINET | organoid | `results/DeepKINET/DeepKINET_organoid.h5ad` | Same layer keys | E, F gamma |
| TSvelo | RPE1 alpha | `results/TSvelo/rpe1/l0_TSvelo.h5ad` and `results/TSvelo/rpe1/l1_TSvelo.h5ad` | `alpha` from both lineages | B alpha, C alpha |
| TSvelo | RPE1 stream | `results/TSvelo/rpe1/TSvelo.h5ad` | `velocity` from the merged result | A |
| Dynamo(label data) | RPE1 | `indicators/gene/scEU_rpe1_analyze/dynamo/dynamo_rpe1.h5ad` | `cell_wise_alpha`; labelled-data stream state | A, B alpha, C alpha |

The selected DeepVelo files must contain `velocity_unspliced`, identifying
the unspliced-prediction result family. DeepKINET's formal seed-123 files
supply gamma only. The merged TSvelo file has no alpha layer and cannot
replace the two lineage files for parameter recovery. Organoid TSvelo and
Dynamo outputs are not added to this experiment.

Reference tables remain `references/kinetics/rpe1_table_s1.csv`,
`references/kinetics/rpe1_alpha_genes.csv`, and
`references/kinetics/organoid_table_s4.csv` under the kinetics artifact root.

## Processing and scoring

### RPE1

Direct-layer RPE1 inputs contain 2,793 pulse-subset cells.

1. Group the configured alpha/gamma layers by `Cell_cycle_possition`, take
   means, and reindex positions 0--300.
2. Intersect genes across the available method parameters and their
   references, preserving reference-table order.
3. Write the grouped raw matrix to CSV and reload it.
4. Only for DeepVelo RPE1, clip negative grouped values to zero.
5. Apply a centered 20-position moving mean with `min_periods=1`.
6. Divide each gene by the median of its own smoothed trajectory; replace a
   zero median with 1. Transform as `log2(ratio + 1e-6)`.
7. Restore the original missing-position mask; align prediction and
   reference genes and timepoints.
8. Write and reload the processed trajectory; apply the existing CSV
   round-trip to the aligned reference before scoring.
9. Compute Spearman correlation on paired finite values. Fewer than two
   valid pairs or constant trajectories yield an undefined score.

Method-specific extraction is retained:
- Dynamo alpha is read from its CSR layer, binned using the current integer
  position boundaries, and stripped of all-zero genes. Its resulting
  timepoint coverage must not be padded silently to mimic direct-layer inputs.
- TSvelo alpha is averaged by cell ID across overlapping lineage cells,
  joined to cycle positions, and then averaged by position and reindexed.
  The position label is retained from the first occurrence of each cell.

Panel B applies a display-only z-score separately to processed reference
and prediction profiles. This z-score is not used for Spearman or turnover
classification. Panel C displays selected-gene Spearman scores.

Panel D independently classifies reference and predicted processed
alpha/gamma trajectories: cosine <= -0.5 is cooperative, cosine >= 0.5 is
destabilizing, and intermediate values are neutral. Genes without finite
paired values or with zero vector norm are unclassified. Only DeepVelo and
cellDancer provide both parameters in this comparison. Save both class-count
and per-gene membership tables.

### Organoid

The source objects contain 3,831 cells before branch-specific selection.

| Branch | Exact Monocle label | SOM clusters | Aligned time |
| --- | --- | --- | --- |
| 1 | `1` | 4, 6, 8, 9, 10, 11 | `3.5 - monocle_pseudotime` |
| 2 | `2` | 1, 2, 6, 7, 8, 11 | `monocle_pseudotime - 6.5` |

Exclude negative aligned times and retain zero. Use `pd.cut` with 200
equal-width bins over the retained aligned-time range and
`include_lowest=True`; average values within bins and reindex all 200 bins.
Smooth with a centered window of 5 and `min_periods=1`. Divide each gene by
the mean of its first 50 smoothed bins plus `1e-6`. Do not clip negative
parameter values or apply log2 or display z-score on this path.

Align genes and timepoints with each branch reference, preserve the
processed-trajectory and aligned-reference CSV round-trips, and calculate
paired-finite Spearman scores. Panel F selects the configured branch genes;
report actual available genes and counts rather than assuming all configured
genes were scored. Current panel data contain six genes in branch 1 and
eight in branch 2.

The evaluator also exports cosine similarity and valid-sample counts.
Component scores are averaged within systems, averaging the two organoid
branches first, and then equally averaged across available systems.
These full scoring summaries are distinct from the selected-gene panels.

## Stream rendering contract

Panels A and E use gene-level velocity and its low-dimensional projection,
not the binned, smoothed, or normalized parameter trajectories above.

For cellDancer, DeepVelo, and merged TSvelo, validate the saved
`velocity_umap` against `X_umap`, cell order, dimensions, and finite values.
For DeepKINET validate `DeepKINET_velocity_umap` and its recorded projection
provenance; its graph keys are `DeepKINET_velocity_graph` and
`DeepKINET_velocity_graph_neg`. Key presence alone is not proof of provenance.
Reuse validated saved projections; do not unconditionally recompute moments,
neighbors, graphs, or embeddings at plotting time. A necessary recomputation
must state its expression layer, neighbor parameters, seed, software versions,
and reason before replacing a panel.

RPE1 PNG export uses `bbox_inches="tight"` and `pad_inches=0.1` at 600 dpi, matching the retained Dynamo panel canvas. This changes only surrounding whitespace, not the saved embedding, velocity projection, or any scores.

DeepKINET RPE1 can obtain labels, palette, and the comparison embedding from
`results/scvelo_deterministic/scv_deterministic_rpe1.h5ad`, aligned by cell ID.
This object does not supply DeepKINET velocity. Verify embedding compatibility
before reusing or recomputing a projection. Organoid labels use the configured
cluster/cell-type annotation. Retain Dynamo's labelled-data rendering route.

## Current figure and numerical scope

The current panels use frozen central-result inputs. C contains the 19 genes
available after evaluator intersection, not all 32 configured genes. F contains
six genes in branch 1 and eight in branch 2; the displayed numbers are means.
These selected-gene statistics must not be replaced by full-gene summaries.
DeepKINET's current plotted Spearman scores agree exactly with the historical
seed-123 bundle for all 19 RPE1 and 6/8 organoid selected genes. Complete-gene
CSV scores differ slightly (cause not established), so seed identity must be established
from H5AD reproducibility metadata and input checksums, not CSV byte identity.

Before replacing any future panel, freeze its inputs, render into a fresh
directory, compare it with an explicitly supplied external archive, and validate
all plot rows and composite content. Assembly-file identity alone does not
establish composite synchronization. DeepKINET's seed-123 designation describes
its own method run, not training every method under one seed.

## Entry points

The evaluator resolves references and `method-results/kinetics/` aliases against
the user-supplied kinetics artifact root. Quantitative plotting resolves paths
against the resulting evaluator output directory; stream rendering resolves
H5AD paths against the velocity root. Each renderer requires a fresh output
directory so a reproduction cannot accidentally overwrite preserved panels.

~~~bash
python evaluation/kinetics_recovery/metrics.py --artifact-root /path/to/kinetics-artifacts --check-only
python evaluation/kinetics_recovery/metrics.py --artifact-root /path/to/kinetics-artifacts --output /path/to/fresh-kinetics-output
python evaluation/kinetics_recovery/plotting/figure5.py --artifact-root /path/to/fresh-kinetics-output --output /path/to/fresh-quantitative-output
python evaluation/kinetics_recovery/plotting/streams.py --artifact-root /path/to/velocity-root --panel rpe1 --output /path/to/fresh-rpe1-streams
python evaluation/kinetics_recovery/plotting/streams.py --artifact-root /path/to/velocity-root --panel organoid --output /path/to/fresh-organoid-streams
python evaluation/kinetics_recovery/plotting/compare_figure5.py --old /path/to/archived-fig5 --new evaluation/kinetics_recovery/figures/fig5 --output /path/to/fresh-comparison
~~~

Use the software versions recorded in the `results/central_results/manifest.json`. The quantitative
renderer also requires `matplotlib-venn` (recorded run: 1.1.2). Scatter jitter
is seeded only for reproducible display; evaluator scores are unaffected.
For another run, supply different fresh output paths rather than rerunning into
these completed directories. The Dynamo A renderer remains unavailable; retain
its recorded manual asset and labelled-data exception until separately resolved.
