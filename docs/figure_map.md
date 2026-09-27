# Figure map

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

Figure 5 is the gene-level real-dataset evaluation. Figure 6 is the
scEU-seq kinetic-parameter recovery evaluation, with active source panels in
`evaluation/kinetics_recovery/figures/fig6/`. Its B–D/F quantitative panels
and A/E streams are mapped in the assembly manifest; E includes the
Dynamo(label data) positive-control stream. The current Figure 6 composite
and its checksums are `figures/figSix.pdf` and
`figures/fig6_composite_manifest.json`. The historical kinetics `fig5`
paths are preserved only under `evaluation/kinetics_recovery/archive/fig5/`.

Figure 7 is the integrated method comparison. Its repository-baseline PDF is
`figures/figSeven.pdf`; the source tables remain under
`evaluation/summary/`, pending a separate Figure 7 audit.
