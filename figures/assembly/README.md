# Paper figure assembly assets

Every materialized panel under this directory is a regular repository file.
`assembly_manifest.tsv` records its paper figure, panel, canonical source,
materialized path, and SHA256; unresolved independent sources remain explicit.

Figure 3 contains three selected gene/cell comparisons. Figure 4 contains
simulation panels; its native schematic/assembly source is retained in
`fig4/authoring/source.pptx`, with the author-export identity verified in its
manifest. A standalone schematic export is not claimed. Its final D
panel is the paired Cell/Gene comparison under `evaluation/cell/simulation/`;
the three gene-only scenario summaries remain supporting analyses. The
approved user-exported composite is recorded in `../fig4_composite_manifest.json`.
Figure 5
contains real-data gene-kinetics panels; the canonical source for its D
overview remains unverified. Figure 6 has 17 materialized source panels,
including Dynamo(label data) and TSvelo in organoid E, with its validated user-exported
composite recorded in `../fig6_composite_manifest.json`. Figure 7's original
editable whole-panel source is `fig7/fig7.pptx`, including ten native charts
and their embedded workbooks. The current author-exported,
cropped `../figSeven.pdf` and its numerical inputs are recorded in
`../fig7_composite_manifest.json`. The old kinetics Figure 5 material is kept
in the archive and has no active assembly aliases.

Supplementary Figure 10's GeoMedian-4 CBDir and BATC source panels are
materialized under `suppfig10/`. Supplementary Figure 11 has two materialized
veloAE source panels. The exact online full-page PDFs for Supplementary Figures 9 and 11 are
retained under the root `figures/` directory. Retaining a final PDF does not
prove independent editable or numerical reconstruction; these boundaries
remain explicit in the composite manifests. Figure 1 also retains the exact
online PDF, while its original authoring source predates the approved three
text corrections. Figure 5D left/right overview images are recovered from the author slide at
the D label position. The formerly misidentified branching-gene summary is
retained only in `figures/archive/fig5/`.

Supplementary Figures 12 and 13 are registered as whole-composite materials.
S12 is the approved simulation Cell composite with its Table 7 provenance;
S13 is the approved TiVelo root-control composite with its editable PPTX and
panel manifest. Their detailed provenance remains beside each composite.
