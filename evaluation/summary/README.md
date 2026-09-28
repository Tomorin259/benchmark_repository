# Evaluation summary

This module combines already computed cell-, gene-, and kinetics-recovery
tables into paper-level summaries. It does not reopen method-result H5AD files
or recalculate a lower-level metric. Figure 7 and other cross-experiment
outputs owned by this module are written under its `results/` and `figures/`
directories.

The historically named `results/figure6_summary/` package contains compact
cross-experiment inputs for the paper's Figure 7. The directory and CSV
filenames retain historical figure numbers; they do not select an older
evaluation. Full final figure assembly remains a documented manual step.

The Real-CBDir columns in the three `fig2_gene_metric_*.csv` summary files
come from the reference-30 `new_cbdir` values in
`evaluation/gene/results/gene_cbdir_reference30_all.csv`. Refresh the 63
dataset-method means and their nine equal-dataset method means and sample
standard deviations with:

```bash
python3 evaluation/summary/refresh_real_cbdir_summary.py --write
```

Run the same command without `--write` to check that the committed summaries
still match the reference-30 results. This refresh does not alter simulation,
kinetic-parameter, or cell-level metrics.

## Kinetic-parameter summaries

The kinetic columns in the same package use
`evaluation/kinetics_recovery/results/central_results/gene_metrics.csv`,
including the corrected organoid Branch 1 reference direction and the
labelled organoid Dynamo alpha control. `results/final_seed123/` is historical
provenance, not an active summary input.

Average finite per-gene Spearman scores within each dataset/branch/parameter
component, average the organoid branch means equally, then average available
system means equally. Report the sample standard deviation across systems;
leave it undefined for a single system. Counts are finite gene scores, not
cells, and unavailable method parameters are not filled with zero. These
full-gene summaries differ from Figure 6F's selected-gene annotations.

Refresh the component, system, method, and WPS-transfer tables with:

```bash
python3 evaluation/summary/refresh_kinetics_summary.py --write
```

The refresh also updates only these four tables' checksums in `artifacts.tsv`.
Without `--write`, the command checks the committed tables and checksums. It verifies
the recomputed means/counts against the three canonical `central_results`
aggregation tables. Unrounded values are used throughout; only the WPS
transfer applies the existing `mean + 1` display transformation. No input
scores, cell-level results, simulation results, or figure assets are changed.
Updating these tables does not regenerate the repository-baseline Figure 7
PDF or the author's PowerPoint assembly.
