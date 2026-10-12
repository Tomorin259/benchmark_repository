# Expanded scEU-seq reference evaluation

These supplementary source results are independent of the selected-gene
Figure 6 panels and the canonical Figure 7 aggregation. They do not replace
`../central_results/` and are not an overall method ranking.

## Gene universe

`gene_scope.csv` records every reference row, original order, candidate order,
source selection, and exclusion reason. RPE1 uses the frozen 250-gene overlap
of Table S1 and the conventional inputs (out of 528 reference genes).
Organoid uses the same 270 genes for both branches: the retained overlap of
the 295 differentiation-associated genes selected in Table S3 and profiled
in Table S4. The 25 unavailable selected genes and six housekeeping controls
are explicitly excluded. This is a broader sensitivity analysis, not the
19 RPE1 or 6+8 lineage examples in Figure 6.

Reference: Battich et al., *Science* (2020),
[doi:10.1126/science.aax3072](https://doi.org/10.1126/science.aax3072),
Tables S1, S3 and S4. Original source hashes and portable external locators
are frozen in `manifest.json`; large H5ADs are registered in `artifacts.tsv`.

## Scores and missing values

Canonical scores are read, restricted to this universe, and combined in
memory with `tsvelo_organoid_gene_metrics.csv`. No duplicate canonical
score table is stored here. TSvelo organoid uses the two lineage `alpha`
layers, not the merged velocity-only H5AD. Uppercase symbols are mapped
back to mouse reference symbols only after original gene-ID validation;
`tsvelo_organoid_gene_mapping.csv` contains the 37 checked mappings.
Overlapping cells are averaged; branch metadata must agree. The merged
input contains 3,831 cells. The shared branch processing selects 1,053 and
287 cells, respectively, and retains 37 scoring rows per branch, of which
20 have finite Spearman scores. The other 17 have constant predictions,
not scores of zero. No gamma score is inferred from a gene-level gamma field.

The shared evaluator settings supply the branch masks, 200 time bins,
five-bin smoothing, early 50-bin baseline and paired-finite Spearman scoring.
Branch 1 reference columns are reversed; Branch 2 is unchanged. The raw and
processed CSV round-trips preserve the audited expanded TSvelo calculation.
RPE1 uses the already-frozen canonical trajectories, without retraining.

`coverage_summary.csv` separates candidate, present, finite, undefined and
unavailable counts. Means/medians exclude non-finite scores; unavailable
genes are never zero-filled. `common_gene_summary.csv` additionally compares
only common finite genes, with exact membership in `common_gene_membership.csv`:

- `core_conventional`: cellDancer/DeepVelo for alpha; also DeepKINET for gamma.
- `all_conventional`: adds TSvelo where alpha scores exist.
- `all_displayed_methods`: additionally includes Dynamo where alpha scores exist.

Common gene counts for alpha are 250/65/42 in RPE1 and 270/20/11 in each
organoid branch. Gamma counts are 249 for RPE1 and 270 per organoid branch.
These are separate sensitivity summaries, not substitutions for each
method's available-gene boxplot mean. Dynamo remains an information-rich
label-data comparator; unequal gene and cell coverage still applies.

## Reproduction

Run from the repository root in the kinetics evaluation environment:

```bash
python -m evaluation.kinetics_recovery.expanded_reference --output NEW_SUMMARY_DIRECTORY
python -m evaluation.kinetics_recovery.expanded_reference --artifact-root ARTIFACT_ROOT --output NEW_RECOMPUTATION_DIRECTORY
python -m evaluation.kinetics_recovery.plotting.expanded_reference --output NEW_FIGURE_DIRECTORY
```

The first command only regenerates summaries from committed scores. The
second validates every external input hash, rebuilds gene scope and
recomputes TSvelo. All output directories must be new or empty. The
rendering command reuses `figure6.py::plot_branch` and `plot_trajectory`.
Three 600-dpi PNGs and one ten-page, 50-genes-per-page PDF are retained at
`../../figures/expanded_reference/`; extra plot-data CSVs in fresh rendering
runs are rebuildable and not duplicate committed score sources.

Status: module source materials are ready. A final supplemental PPT/composite,
manuscript figure number and online insertion are deliberately not registered.
