# Paper figure assembly assets

Every materialized panel under this directory is a regular repository file.
`assembly_manifest.tsv` records its paper figure, panel, canonical source,
materialized path, and SHA256; unresolved independent sources remain explicit.

Figure 3 contains three selected gene/cell comparisons. Figure 4 contains
simulation panels; its independent schematic source is missing. Its final D
panel is the paired Cell/Gene comparison under `evaluation/cell/simulation/`;
the three gene-only scenario summaries remain supporting analyses. The
approved user-exported composite is recorded in `../fig4_composite_manifest.json`.
Figure 5
contains real-data gene-kinetics panels; the canonical source for its D
overview remains unverified. Figure 6 has 16 materialized source panels,
including Dynamo(label data) in organoid E, with its validated user-exported
composite recorded in `../fig6_composite_manifest.json`. Figure 7 retains the
repository-baseline `../figSeven.pdf`; no independent whole-panel source is
recorded. The old kinetics Figure 5 material is kept in the archive and has
no active assembly aliases.

Supplementary Figure 10's GeoMedian-4 CBDir and BATC source panels are
materialized under `suppfig10/`. Supplementary Figure 11 has two materialized
veloAE source panels. Its manual full-page composite and Supplementary Figure 9
remain unresolved in the manifest.

Supplementary Figures 12 and 13 are registered as whole-composite materials.
S12 is the approved simulation Cell composite with its Table 7 provenance;
S13 is the approved TiVelo root-control composite with its editable PPTX and
panel manifest. Their detailed provenance remains beside each composite.
