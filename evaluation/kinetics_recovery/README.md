# Kinetic-parameter recovery

Compare saved method outputs with published scEU-seq rate references for RPE1
and intestinal organoid. The evaluator reads H5ADs and reference tables;
method training is a separate step.

## Results and scope

| Material | Purpose |
| --- | --- |
| [results/central_results/](results/central_results/) | Canonical per-gene scores and summaries on the original evaluation scope |
| [figures/fig6/plot_data/](figures/fig6/plot_data/) | Exact selected rows plotted in Figure 6 |
| [figures/fig6/](figures/fig6/) | Current parameter panels and saved velocity fields |
| [results/expanded_reference/](results/expanded_reference/README.md) | Expanded reference scores for Figure 7 and Supplementary Figure 14 |
| [results/joint_recovery/](results/joint_recovery/) | Common-gene alpha and gamma comparison for cellDancer and DeepVelo |
| [results/turnover/](results/turnover/) | Full-precision class F1, confidence intervals and paired tests for Supplementary Table S11 |

Figure 6B/C displays 19 RPE1 genes surviving evaluator intersection from the
configured 32. Figure 6F displays six branch-1 and eight branch-2 genes.
Panel D classifies 250 matched genes with both parameters, separately from the
19-gene display. A/E are cell-level velocity fields.
[Figure 7](../summary/README.md) uses each method’s finite scores within the
expanded reference sets of Supplementary Figure 14 (250 RPE1 genes and 270
organoid genes per branch). Both organoid branches and both biological systems
are equally weighted at their respective aggregation steps. The expanded
TSvelo organoid alpha scores come from the lineage-specific outputs.
Per-method coverage and the separate common-finite-gene comparisons are
available in the expanded-reference directory.

## Joint alpha and gamma recovery

The comparison starts from the public
[per-gene Spearman scores](results/central_results/gene_metrics.csv). Within
RPE1 and each organoid branch, a gene enters the comparison when both
cellDancer and DeepVelo have finite scores for both parameters. This gives
249 RPE1 genes and 270 genes per organoid branch. Each row in
[joint scores](results/joint_recovery/common_gene_scores.csv) contains both
parameter scores for one method and one gene in a system or branch, for
1,578 rows total.

The [comparison summary](results/joint_recovery/common_gene_comparison.json)
records the input SHA-256, common gene counts, component means, system means,
and the descriptive Pareto front for each system. For each method and
parameter, genes contribute equally to the component mean. The organoid
system mean gives its two branch means equal weight. A method dominates when
it has at least as high a mean Spearman correlation on both parameters and a
strictly higher mean on at least one. DeepVelo has higher mean alpha and
gamma correlations for both RPE1 and organoid.

The expanded-reference analysis defines its common finite gene sets
separately for each parameter and each configured method group. The joint
comparison uses one gene set per component for all four method–parameter
pairs.

Run from the repository root with Python 3.8 or later and a new output directory.

```bash
python -m evaluation.kinetics_recovery.joint_recovery --output outputs/kinetics-joint-recovery
```

## RPE1 turnover classification (Supplementary Table S11)

The [full-precision table](results/turnover/class_f1.csv) compares cellDancer
and DeepVelo on the 250 matched RPE1 genes classified in Figure 6D. The script
reads each method's gene memberships and class counts, checks the counts and
the shared scEU-seq reference labels, then computes class F1 and the unweighted
Macro-F1. The difference column is DeepVelo minus cellDancer.

The 95% confidence intervals are the 2.5th and 97.5th percentiles from 10,000
bootstrap resamples of genes within their scEU-seq class. Both methods use the
same resamples. The two-sided test independently swaps the two method labels
within each gene over 50,000 permutations. Holm adjustment covers the three
class tests. Macro-F1 has a confidence interval and blank P-value fields.
The bootstrap and permutation use independent streams from seed `20261111`.

Run from a Git checkout with NumPy (release computation: 1.23.5) and a new
output directory.

```bash
python -m evaluation.kinetics_recovery.turnover --output outputs/kinetics-turnover
```

For the standalone replay, pass the directory containing the four
`D_rpe1_*_kinetic_class*.csv` files from Figure 6D. The repository filenames
without the `D_` prefix also work. The command writes `class_f1.csv` and a
[manifest with input and output checksums](results/turnover/manifest.json).

```bash
python turnover.py --input-dir ../../main_figures/Figure_6 --output s11-reproduction
```

## Evaluate and render

Run from the repository root. Obtain the [original reference tables](../../data/UPSTREAM_REFERENCES.md)
and [saved method outputs](../../data/RESULT_DOWNLOADS.md), with software
versions from [the current run manifest](results/central_results/manifest.json).
The quantitative renderer also requires `matplotlib-venn` (recorded run: 1.1.2).
Use the [kinetics dependency record](../../docs/environments.md#evaluation)
for numerical evaluation and the plotting dependencies for figures.

`--artifact-root` is the local directory containing both
`references/kinetics/` and `method-results/kinetics/`. The
[run manifest](results/central_results/manifest.json) pairs each method-result
alias with its original velocity-root path and SHA-256; expose the downloaded
file at the alias consumed by `config.json`. Reference filenames are listed
below. [Local roots](../../data/README.md#local-data-roots) distinguishes this
layout from the Cell and Gene roots.

```bash
python evaluation/kinetics_recovery/metrics.py \
  --artifact-root /path/to/kinetics-artifacts --check-only
python evaluation/kinetics_recovery/metrics.py \
  --artifact-root /path/to/kinetics-artifacts --output outputs/kinetics-evaluation
python evaluation/kinetics_recovery/plotting/figure6.py \
  --artifact-root outputs/kinetics-evaluation --output outputs/kinetics-panels
python evaluation/kinetics_recovery/plotting/streams.py \
  --artifact-root /path/to/velocity-root --panel rpe1 --output outputs/rpe1-streams
python evaluation/kinetics_recovery/plotting/streams.py \
  --artifact-root /path/to/velocity-root --panel organoid --output outputs/organoid-streams
```

Use fresh output directories. The evaluator resolves reference and
`method-results/kinetics/` aliases under the kinetics artifact root.
`--check-only` reads the H5ADs and reference tables, computes the component
scores and prints them without saving. `--methods` selects a method subset;
the default evaluates every configured method.

Numerical evaluation writes `gene_metrics.csv`, `component_means.csv`,
`system_means.csv`, `method_summary.csv` and processed matrices under
`matrices/`. Check the method/parameter/system/branch fields and finite-gene
counts when comparing with `results/central_results/`.
Quantitative plotting reads evaluator output; streams read saved velocity-root
H5ADs. [Plotting details](plotting/README.md) describe display settings and
[expanded rendering](results/expanded_reference/README.md).

The Figure 6 renderer's input paths are in
`plotting/figure6_config.json`. The default complete evaluation supplies its
configured matrices and scores; a method-subset evaluation supplies only that
subset. `streams.py` takes the velocity root separately and reads the stored
velocity fields shown in the native-input table.

## Selected native inputs

Paths are relative to the velocity root; configuration aliases may expose the
same identities under the kinetics artifact root.

| Method / task | Saved H5AD | Consumed fields |
| --- | --- | --- |
| cellDancer RPE1 / organoid | `results/celldancer/celldancer_{rpe1,organoid}.h5ad` | `celldancer_alpha`, `celldancer_gamma`, `velocity` |
| DeepVelo RPE1 / organoid | `results/deepvelo/results_with_vu/deepvelo_{rpe1,organoid}.h5ad` | `cell_specific_alpha`, `cell_specific_gamma`, `velocity`; `velocity_unspliced` identifies the selected family |
| DeepKINET RPE1 / organoid | `results/DeepKINET/DeepKINET_{rpe1,organoid}.h5ad` | `degradation_rate` (gamma), `DeepKINET_velocity` |
| TSvelo RPE1 alpha | `results/TSvelo/rpe1/l0_TSvelo.h5ad`, `l1_TSvelo.h5ad` | Lineage `alpha` |
| TSvelo RPE1 / organoid streams | `results/TSvelo/rpe1/TSvelo.h5ad`, `results/TSvelo/organoid_ori/TSvelo.h5ad` | Merged `velocity` and saved projection; no selected-gene parameter scores from the merged files |
| Dynamo labelled RPE1 | `indicators/gene/scEU_rpe1_analyze/dynamo/dynamo_rpe1.h5ad` | `cell_wise_alpha` and labelled-data stream |
| Dynamo labelled organoid | `indicators/gene/scEU_organoid_analyze/dynamo/pulse120_one_shot_corrected_20260926/dynamo_pulse120_corrected.h5ad` | Alpha-only 120-min pulse control; alias `method-results/kinetics/dynamo/organoid_labelled.h5ad` |

Dynamo's labelled inputs are distinct from the label-free RNA-seq configuration.
Organoid Dynamo uses 1,373 labelled cells over a 3,831-cell background; the
3,831-cell, 153-gene merged TSvelo stream reuses the reference UMAP.
Reference files under the kinetics root are
`references/kinetics/rpe1_table_s1.csv`, `rpe1_alpha_genes.csv` and
`organoid_table_s4.csv`.

## Processing and scoring

### RPE1

Direct-layer inputs contain 2,793 pulse-subset cells.

1. Average by `Cell_cycle_possition`, reindex positions 0–300, and intersect
   prediction/reference genes in reference order.
2. Write/reload the grouped CSV; clip negative grouped values only for DeepVelo.
3. Apply a centered 20-position moving mean (`min_periods=1`).
4. Divide each gene by its smoothed median (use 1 for zero), then take
   `log2(ratio + 1e-6)`.
5. Restore the missing-position mask, align to the reference and preserve the
   processed/reference CSV round-trips before paired-finite Spearman scoring.
   Fewer than two pairs or constant profiles yield undefined scores.

Dynamo alpha uses CSR extraction and integer position bins, drops all-zero
genes and retains its actual timepoint coverage. TSvelo alpha averages
overlapping lineage cells by cell ID before position averaging/reindexing;
the first recorded position label is retained.

B z-scores prediction/reference curves independently for display only.
C reports finite displayed-gene means and counts.
D independently classifies processed reference and predicted alpha/gamma:
cosine ≤ -0.5 cooperative, ≥ 0.5 destabilizing, otherwise neutral.
Missing pairs and zero norms are unclassified. Only cellDancer and DeepVelo
supply both parameters. Class F1 is `2 TP / (2 TP + FP + FN)`;
Macro-F1 is the unweighted three-class mean on the same 250-gene universe.
Counts and gene memberships are exported; the displayed values are
0.324 (cellDancer) and 0.340 (DeepVelo).

### Organoid

Select branches from the 3,831-cell source:

| Branch | Monocle label | SOM clusters | Aligned time |
| --- | --- | --- | --- |
| 1 | `1` | 4, 6, 8, 9, 10, 11 | `3.5 - monocle_pseudotime` |
| 2 | `2` | 1, 2, 6, 7, 8, 11 | `monocle_pseudotime - 6.5` |

Exclude negative aligned times and retain zero. Use `pd.cut` with 200
equal-width bins and `include_lowest=True`, average/reindex all bins,
and smooth with a centered window of 5 (`min_periods=1`).
Divide by each gene's first-50-smoothed-bin mean plus `1e-6`.
This path does not clip negative rates, log-transform or display-z-score them.

The Table S4 export's branch-1 columns run opposite to the evaluator convention:
reverse only those reference columns before assigning timepoints 0–199.
Branch-2 references and prediction trajectories are unchanged.
Align genes/timepoints, retain the CSV round-trips, then score paired finite
values. Export cosine similarity and valid-sample counts as well.

### Aggregation

Selected-gene means use available finite scores within each displayed branch;
an organoid display average weights the two branch means equally.
Full summaries average finite gene scores within components, average organoid
branches equally, then average available system means equally.
Missing parameters are not zero-filled. These statistics have different gene
sets and must remain separate.

The DeepKINET H5AD uses seed 123.

## Saved velocity fields and final assembly

A/E use saved cell-level velocity projections, separately from processed rate
trajectories. Stream rendering checks cell IDs, finite values, shapes and the
configured comparison embedding. DeepKINET checks also include recorded graph
and projection provenance.

DeepKINET RPE1 may obtain aligned labels/palette from
`results/scvelo_deterministic/scv_deterministic_rpe1.h5ad`; DeepKINET velocities
come from its own saved H5AD. Rendering manifests record input identities and
validation checks; records of the original projection computation are retained
where available.
`streams.py --panel rpe1` renders all five Figure 6A methods with the same
scVelo call on the shared RFP/GFP layout. The configured Dynamo fields are
mapped to the standard plotting names in memory. Labelled organoid Dynamo
retains its original image asset. [Source Data](../../data/SOURCE_DATA.md)
includes the saved coordinates and velocities.

The final PDF and panel bindings are in the [figure map](../../docs/figure_map.md).
