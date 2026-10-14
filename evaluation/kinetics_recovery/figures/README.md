# Kinetics figure assets

`fig6/` is the active Figure 6 panel directory for scEU-seq kinetic-parameter recovery. B–D/F PNGs and their plotted rows are directly under `fig6/` and `fig6/plot_data/`; RPE1 A and organoid E streams are under `fig6/figure6A_rpe1/` and `fig6/figure6E_organoid/`. Dynamo(label data) is a manually assembled positive-control stream in both A and E, with source provenance recorded beside E.

Quantitative panels use `../results/central_results/`. The assembly copies and hashes are under `../../../figures/assembly/fig6/` and its manifest. The current composite is `../../../figures/figSix.pdf`; `../../../figures/fig6_composite_manifest.json` records the source-panel and composite checksums.

`../archive/fig5/` retains unique historical panels and provenance. The complete
former kinetics Figure 5 bundle is available at
[Git revision `bed37b629dbc966109c8cb92dd22aca0810caf2f`](https://github.com/Tomorin259/benchmark_repository/tree/bed37b629dbc966109c8cb92dd22aca0810caf2f/evaluation/kinetics_recovery/archive/fig5);
see the owning module README for comparison materialization. Current manuscript
Figure 5 is a separate gene-level figure.

`expanded_reference/` contains three reviewed 600-dpi boxplot PNGs and one
10-page RPE1 250-gene trajectory PDF. They are supplementary source materials,
not replacements for `fig6/` and not a registered paper composite. Their
numerical scope and hashes are owned by `../results/expanded_reference/manifest.json`.
