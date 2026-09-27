# Paper figure assembly assets

Every materialized panel under this directory is a regular repository file.
`assembly_manifest.tsv` records its paper figure, panel, canonical source,
materialized path, and SHA256; unresolved independent sources remain explicit.

Figure 3 contains three selected gene/cell comparisons. Figure 4 contains
simulation panels; its independent schematic source is missing. Figure 5
contains real-data gene-kinetics panels; the canonical source for its D
overview remains unverified. Figure 6 has 16 materialized source panels,
including Dynamo(label data) in organoid E, with its online-matched composite
recorded in `../fig6_composite_manifest.json`. Figure 7 retains the
repository-baseline `../figSeven.pdf`; no independent whole-panel source is
recorded. The old kinetics Figure 5 material is kept in the archive and has
no active assembly aliases.

Supplementary Figure 10's GeoMedian-4 CBDir and BATC source panels are materialized under `suppfig10/`. Supplementary Figures 9 and 11 are recorded as unresolved final-composite sources in the manifest.
