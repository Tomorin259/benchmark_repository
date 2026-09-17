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

Figure 5 uses the sole canonical panel directory
`evaluation/kinetics_recovery/figures/fig5/` and the current central-result
evaluation. Its composite mapping and JPEG-content validation are in
`figures/fig5_composite_manifest.json`. Archived figures are external assets.
