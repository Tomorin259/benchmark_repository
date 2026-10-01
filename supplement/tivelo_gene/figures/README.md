# Final TIvelo materials for Supplementary Fig. 13

See [the module README](../README.md) for inputs and rendering commands.
This directory owns the approved A/B compact material, simulation C and bone
marrow D, source panels, plotted data and provenance. The final paper composite
and editable PPT belong in `figures/assembly/suppfig13/` at repository root.

A/B rows are branching and hybrid. Columns are root selection, known order,
main-path DPT and secondary-path DPT. Known bins are divided by nine for a
relative 0--1 progression scale. Each path's DPT is independently normalized;
the common colorbar does not imply a shared elapsed-time unit. Gray cells are
outside the respective DPT path. Coordinates, colors and root markers are
unchanged. The compact assembler crops the original gradient from
`ab_assets/order_comparison_material.pdf`.

D uses native full-data `obs.dpt_pseudotime` from the prespecified first paired
seed, 20261001, on identical archived t-SNE coordinates for all 5,780 cells.
It is not Palantir-colored or newly computed UMAP. Both maps use a fixed 0--1
color range. The metrics show all three paired optimizer seeds. Path-wise
orientation DPT in B and full-data DTI DPT in D are different stages. The
real-data intervention changes the child graph and directed neighbors as well
as the root.

Historical source manifests describe the original renderers. Consolidated
scripts share identical styling through `scripts/plot_style.py`; no fitting,
score calculation, source data or approved panel binaries changed in cleanup.
