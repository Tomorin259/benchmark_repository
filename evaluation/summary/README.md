# Evaluation summary

This module combines already computed cell-, gene-, and kinetics-recovery
tables into paper-level summaries. It does not reopen method-result H5AD files
or recalculate a lower-level metric. Figure 7 and other cross-experiment
outputs owned by this module are written under its `results/` and `figures/`
directories.

The `results/figure7_summary/` package contains the canonical compact Figure 7
inputs. Its directory, CSV names and method-label column use Figure 7 names.
Final assembly uses the original editable PowerPoint and native chart workbooks;
PDF export remains a documented manual step.

The Real-CBDir columns in the three `fig7_gene_metric_*.csv` summary files
come from the reference-30 `new_cbdir` values in
`evaluation/gene/results/gene_cbdir_reference30_all.csv`. Refresh the 63
dataset-method means and their nine equal-dataset method means and sample
standard deviations with:

```bash
python3 evaluation/summary/refresh_real_cbdir_summary.py --write
```

Run the same command without `--write` to check that the committed summaries
still match the reference-30 real results. This refresh does not alter simulation,
kinetic-parameter, or Cell-level metrics. The two simulation Gene columns retain
the frozen `evaluation/gene/results/simulation/gene_cbdir_cossim_reference.csv`
values, summarized within each scenario and then equally across scenarios.

Read-only validation compares unrounded Real-CBDir means and standard
deviations with absolute tolerance `1e-12`, matching the source consistency
check. This permits last-bit floating-point differences between Python
versions, not displayed-score rounding. Coverage, counts, labels, other
metrics, and four-decimal WPS-transfer fields remain exact checks. Validation
does not rewrite the CSVs.

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
Updating these tables does not regenerate the Figure 7 PDF or native chart
workbooks automatically.

## Final Figure 7 assembly

The desktop Gene/kinetics workbook and native Figure 7 charts were synchronized
on 2026-10-02; the Cell workbook already matched the canonical Cell table.
`figures/assembly/fig7/fig7.pptx` contains the original layout, ten native editable
charts and their embedded source workbooks. `figures/figSeven.pdf` is the author's
2026-10-02 WPS re-export, cropped to remove the lower blank slide area.
`figures/fig7_composite_manifest.json` records the inputs, refresh scripts,
assembly, export and PDF hashes and the verification scope.

To reproduce the final export, open the PowerPoint in WPS Office and export
its single slide as a vector PDF. Retain chart sizes, positions, colors and axes.
Set both MediaBox and CropBox to `[0, 315.395, 540, 779.8]` points. Native chart
Excel sheet/range names retain their original bindings; repository CSV names
use Figure 7. This manual export is the final assembly step, not an independent
Python renderer.
