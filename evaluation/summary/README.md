# Evaluation summary

This module combines already computed cell-, gene-, and kinetics-recovery
tables into paper-level summaries. It does not reopen method-result H5AD files
or recalculate a lower-level metric. Figure 6 and other cross-experiment
outputs owned by this module are written under its `results/` and `figures/`
directories.

The historically named `results/figure6_summary/` package contains compact
cross-experiment inputs for the paper's summary figure, including the seed-123
DeepKINET gamma summaries. Full final figure assembly remains a documented
manual step.

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
