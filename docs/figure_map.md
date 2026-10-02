# Figure map

Supplementary Table S6 (veloVI gene-wise coordinate sensitivity) is backed by
the full-precision CSV, transformation script, and manifest in
[supplement/velovi/](../supplement/velovi/). Its A/C post-hoc diagnostics are
separate from the primary cell-level benchmark (condition B).

Scientific panels are generated beside the experiment that owns their data and
metric definition:

- cell-level panels: [`../evaluation/cell/`](../evaluation/cell/);
- gene-level real and simulation panels: [`../evaluation/gene/`](../evaluation/gene/);
- kinetics-recovery panels: [`../evaluation/kinetics_recovery/`](../evaluation/kinetics_recovery/);
- integrated summaries: [`../evaluation/summary/`](../evaluation/summary/);
- ensemble-specific panels: [`../ensemble/`](../ensemble/).

The owning module documents its renderers, inputs, and manual-asset boundaries.
Paper-level composite PDFs are stored in [`../figures/`](../figures/).

[`../figures/assembly/assembly_manifest.tsv`](../figures/assembly/assembly_manifest.tsv)
is the authoritative panel-to-source map. It records panel identity, status,
source and materialized paths, SHA256, and unresolved source gaps. Do not copy
its inventory into prose; derive current panel counts and missing-source lists
from the manifest.

Figure 3 compares selected gene-level and cell-level transition cases; its
three source panels are under `evaluation/gene/figures/panels/selected_edges/`.
Figure 4 is the simulation evaluation, with its final paired Cell/Gene D source
at `evaluation/cell/simulation/figures/simulation_cross_level_cbdir.png`.
Gene-only scenario summaries are supporting analyses, not the final D panel.
Figure 5 is the real-dataset
gene-level evaluation. Their active PDFs are `figures/figThree.pdf`,
`figures/figFour.pdf`, and `figures/figFive.pdf`. Figure 4 is the approved
user-exported, cropped composite; `figures/fig4_composite_manifest.json`
records its checksum and the verification boundary. Its online byte identity
was not independently rechecked during this synchronization.

Figure 6 is the
scEU-seq kinetic-parameter recovery evaluation, with active source panels in
`evaluation/kinetics_recovery/figures/fig6/`. Its B–D/F quantitative panels
and A/E streams are mapped in the assembly manifest; E includes the
Dynamo(label data) positive-control stream. The current Figure 6 composite
and its checksums are `figures/figSix.pdf` and
`figures/fig6_composite_manifest.json`. The historical kinetics `fig5`
paths are preserved only under `evaluation/kinetics_recovery/archive/fig5/`.
Method-colored Figure 6F boxplots inherit the benchmark palette from
`evaluation/cell/plots.py`. Figure 6B trajectory lines retain a documented
high-contrast exception; Figure 6A/E stream colors encode cell clusters.

Figure 7 is the integrated method comparison. Its current author-exported,
cropped PDF is `figures/figSeven.pdf`; canonical summary tables are under
`evaluation/summary/results/figure7_summary/`. The original editable assembly
is `figures/assembly/fig7/fig7.pptx`, including ten native charts and embedded
workbooks. `figures/fig7_composite_manifest.json` links numerical inputs, summary
scripts, assembly, export/crop and output hashes. Final export remains a manual
WPS step; the verification scope is recorded in that manifest.

The online manuscript assigns the matched cell/gene comparison to
Supplementary Figure 9 and GeoMedian-4 evaluation to Supplementary Figure 10.
Figure 10 source panels are in `figures/assembly/suppfig10/`; Figure 11
veloAE A/B source panels are in `figures/assembly/suppfig11/`. Supplementary
Figure 12 is the approved simulation Cell composite at
`figures/suppFigTwelve.pdf`, with numerical provenance in
`figures/suppfig12_composite_manifest.json`. Supplementary Figure 13 is the
approved TiVelo root-control composite under
`figures/assembly/suppfig13/`, with its PDF, editable PPTX, and assembly
manifest. Both are registered as whole composites in the root assembly
manifest.
