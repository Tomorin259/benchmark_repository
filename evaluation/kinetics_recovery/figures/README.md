# Kinetics figure assets

`fig5/` is the sole current Figure 5 panel directory. B--D/F PNGs and exact
plot rows are stored directly under `fig5/` and `fig5/plot_data/`; A/E streams
are in `figure5A_rpe1/` and `figure5E_organoid/`, with source/output manifests.
Dynamo A is the retained labelled-data manual exception.

All quantitative panels consume `../results/central_results/`. DeepKINET is
seed 123 and supplies gamma only. C displays 19 available genes; F displays
six/eight branch genes, with mean annotations. Full-gene summaries have a
different statistical scope. See [the module contract](../README.md).

Old panels and comparison images are archived externally. Current assembly
copies are under `../../../figures/assembly/fig5/`, with canonical source paths
and SHA256 in the assembly manifest. `../../../figures/figFive.pdf` is the
reviewed current composite; `../../../figures/fig5_composite_manifest.json`
records image-content validation of its JPEG-encoded panel images.
