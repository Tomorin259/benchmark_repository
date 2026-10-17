# Figure map

Supplementary Figure S15 and Table S9 are backed by
[`supplement/multivelo/`](../supplement/multivelo/). S15A uses three phase
portraits and a short progression key; S15B uses the single paired-score row
for the two branching subsets. The current author-exported composite and its
editable assembly are registered in
[`figures/suppfig15_composite_manifest.json`](../figures/suppfig15_composite_manifest.json).
The compact four-row S9 LaTeX source is generated from the frozen full-precision
condition summary, not transcribed from the plotted points. The full MURK,
smoothing, complete-trajectory and paired-decrease records remain in the module.

Supplementary Table S6 (veloVI gene-wise coordinate sensitivity) is backed by
the full-precision CSV, transformation script, and manifest in
[supplement/velovi/](../supplement/velovi/). Its A/C post-hoc diagnostics are
separate from the primary cell-level benchmark (condition B).

Figure 1 and Supplementary Figures 9/11 retain the exact online PDF bytes in
`figures/`, checked by online git-blob identity. Authoring candidates and their
verification boundaries are under `figures/assembly/*/authoring/`. Current
Figure 4 has a verified author-export lineage; its standalone A export,
Figure 1 corrected editable A/B must not be inferred from the presence of a
PPTX. Figure 5D sources are identified by native image placement; its old
misidentified gene summary is archived without changing the current PDF.

`figures/export_materials.py` validates sources and generates fresh material
bundles from the map; duplicate source-panel copies are not tracked.

Scientific panels are generated beside the experiment that owns their data and
metric definition:

- cell-level panels: [`../evaluation/cell/`](../evaluation/cell/);
- gene-level real and simulation panels: [`../evaluation/gene/`](../evaluation/gene/);
- kinetics-recovery panels: [`../evaluation/kinetics_recovery/`](../evaluation/kinetics_recovery/);
- integrated summaries: [`../evaluation/summary/`](../evaluation/summary/);
- ensemble-specific panels: [`../ensemble/`](../ensemble/).

The owning module documents its renderers, inputs, and manual-asset boundaries.
Paper-level composite PDFs are stored in [`../figures/`](../figures/).

Gene comparison panels and score summaries use the repository-owned
`evaluation/gene/plot_gene_assets.py` entry point. It consumes the registered
H5AD inputs and committed CSVs; `--summaries-only` consumes CSVs alone. Export
into a fresh directory under `outputs/` before any separately reviewed asset
promotion. It does not rebuild retained raw phase plates or final manual PDF
assemblies. See the Gene README for the exact layer and score-source contract.

[`../figures/assembly/assembly_manifest.tsv`](../figures/assembly/assembly_manifest.tsv)
is the authoritative panel-to-source map. It records panel identity, status,
source and on-demand export paths, SHA256, and unresolved source gaps. Do not copy
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
`figures/fig6_composite_manifest.json`. Unique historical kinetics `fig5`
panels and provenance remain under `evaluation/kinetics_recovery/archive/fig5/`.
The complete historical bundle, including exact copies omitted from the working
tree, is available at
[Git revision `bed37b629dbc966109c8cb92dd22aca0810caf2f`](https://github.com/Tomorin259/benchmark_repository/tree/bed37b629dbc966109c8cb92dd22aca0810caf2f/evaluation/kinetics_recovery/archive/fig5).
Method-colored Figure 6F boxplots inherit the benchmark palette from
`evaluation/cell/plots.py`. Figure 6B trajectory lines retain a documented
high-contrast exception; Figure 6A/E stream colors encode cell clusters.

Figure 7 is the integrated method comparison. Its verified author-exported,
cropped PDF is `figures/figSeven.pdf`, including the corrected Cell CBDir.
Canonical summary tables are under
`evaluation/summary/results/figure7_summary/`. The original editable assembly
is `figures/assembly/fig7/fig7.pptx`, including ten native charts and embedded
workbooks. `figures/fig7_composite_manifest.json` links numerical inputs, summary
scripts, assembly, export/crop and output hashes. Final export remains a manual
WPS step; the verification scope is recorded in that manifest.

Figures 2, 7 and S10 have verified cropped author exports and updated native
assemblies. Figure 2C's CBDir panel and Figure 7's
Cell CBDir column use the corrected Forebrain 30-PC scores; S10 uses the updated
Cell distribution and shared-frame GeoMedian-4 BATC. Figure 2A/B and the other
four Cell metric panels are unchanged. Their composite manifests bind the
current PDF files to the source panels and author exports.

The online manuscript assigns the matched cell/gene comparison to
Supplementary Figure 9 and GeoMedian-4 evaluation to Supplementary Figure 10.
Supplementary Figure 10 SVG sources are in `figures/assembly/suppfig10/`;
Figure 11 veloAE A/B sources are owned by `supplement/veloae/`. Their export
locations are recorded in the root assembly manifest. Supplementary
Figure 12 is the approved simulation Cell composite at
`figures/suppFigTwelve.pdf`, with numerical provenance in
`figures/suppfig12_composite_manifest.json`. Supplementary Figure 13 is the
approved TiVelo root-control composite under
`figures/assembly/suppfig13/`, with its PDF, editable PPTX, and assembly
manifest. Both are registered as whole composites in the root assembly
manifest.
