# Figure 7 summaries

[results/figure7_summary/](results/figure7_summary/) contains compact inputs for
the integrated comparison. These scripts aggregate existing scores without
opening method H5ADs or recalculating metrics.

## Check or refresh tables

Run from the repository root. The default mode checks committed tables
without writing:

```bash
python3 evaluation/summary/refresh_cell_summary.py
python3 evaluation/summary/refresh_real_cbdir_summary.py
python3 evaluation/summary/refresh_sim_cbdir_summary.py
python3 evaluation/summary/refresh_kinetics_summary.py
```

Add `--write` only for an intended refresh. Kinetics refresh also updates the
four affected table checksums in `artifacts.tsv`.
Neither mode regenerates native chart workbooks or the final PDF.

These checks use the standard library and read the committed evaluation
tables listed below. Cell's refresher accepts `--metrics` and `--output`;
the other refreshers declare their source and destination paths in the
script. Select a different evaluator output directory there before
summarizing a new result bundle.
Each refresher expects the full registered comparison scope for its summary.

Successful checks compare the derived rows, coverage and formatted transfer
values with the saved summary tables. Written Figure 7 tables live under
`results/figure7_summary/`. Use their source and count fields to identify the
aggregation represented by each row.

## Sources and aggregation

| Summary | Source | Aggregation |
| --- | --- | --- |
| Five Cell metrics | [280 Cell rows](../cell/results/cell_metrics.csv) | Equal weights for 14 datasets per method |
| Real Gene CBDir | [reference-30 Gene table](../gene/results/gene_cbdir_reference30_all.csv), `new_cbdir` | Available-gene means per dataset, then equal weights for seven datasets |
| Simulated Gene CBDir | Same `new_cbdir` table | Means for 10/20/40 special genes in MURK/branching/hybrid, then equal weights for three scenarios |
| Simulation CosSim | [CosSim table with restored veloVI units](../gene/results/simulation/gene_cbdir_cossim_reference.csv) | Same scenario weighting; historical CBDir in this table is not used |
| Kinetic parameters | [canonical per-gene scores](../kinetics_recovery/results/central_results/gene_metrics.csv) | Finite-gene means per component, equal organoid-branch means, then equal available-system means |

Cell SD uses `ddof=1`; quartiles use linear interpolation
(Excel `QUARTILE.INC`). Real Gene summaries retain sample SD across dataset
means. Kinetics reports sample SD across systems, undefined for a single system.
Counts are finite gene scores, not cells; unavailable parameters are not zero.
Kinetics includes the corrected organoid reference direction and labelled
Dynamo alpha, not historical `final_seed123` or the selected Figure 6 genes.

Forebrain Cell CBDir uses all 30 shared PCs; BATC uses shared UMAP2.
Real/Sim Gene validation permits absolute tolerance `1e-12` for unrounded
means and SDs. Counts, labels, coverage and four-decimal transfer fields
remain exact checks.
CBDir and CosSim are aligned by dataset/method/gene across 630 rows.
Aggregation uses the stored source values before display formatting.

## Final chart assembly

[fig7.pptx](../../figures/assembly/fig7/fig7.pptx) retains ten editable charts and
their embedded workbooks. Native charts apply the existing `mean + 1`
display transform, with unshifted SDs; statistical summaries use the original
scale. The website reads unshifted values.

To reproduce [figSeven.pdf](../../figures/figSeven.pdf), export the single slide
as a vector PDF in WPS Office, retaining positions, colors and axes. Set
MediaBox and CropBox to `[0, 315.395, 540, 779.8]` points.
Workbook sheet/range bindings retain their original names.
[The composite manifest](../../figures/fig7_composite_manifest.json) identifies
current source tables, scripts, author assembly and export.
