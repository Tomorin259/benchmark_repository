# Evaluation summary

This module combines already computed cell-, gene-, and kinetics-recovery
tables into paper-level summaries. It does not reopen method-result H5AD files
or recalculate a lower-level metric. Figure 6 and other cross-experiment
outputs owned by this module are written under its `results/` and `figures/`
directories.

The current `results/figure6_summary/` package contains the refreshed compact
Figure 6 inputs, including the seed-123 DeepKINET gamma summaries. A future
summary plotting module will consume these inputs; full final figure assembly
remains a documented manual step.
