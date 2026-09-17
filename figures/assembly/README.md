# Paper figure assembly assets

Every materialized panel under this directory is a regular repository file; no
panel is a symlink and no panel depends on a path outside the repository.
`assembly_manifest.tsv` records the paper panel, repository-relative canonical
source, materialized path, SHA256, and any unresolved provenance note.

Current assembly: 110 materialized panels, 1 panel with an unverified missing
canonical source (Fig.4D), and 4 missing independent panel sources (Fig.1A,
Fig.1B, Fig.3A, and Fig.6 whole). SupFig.09 has materialized
RowGeometricMedian-4 CBDir and BATC metric panels. These unresolved items are recorded explicitly;
no replacement image was inferred. The preserved top-level PDFs
`../figTwo.pdf` and `../figSix.pdf` remain unchanged; `../figFive.pdf` uses the current central-result panels.

Figure 5 materialized files and hashes match the current canonical panels.
The separately stored composite-content validation is
`../fig5_composite_manifest.json`; matching assembly hashes alone is not proof
of composite contents. The [kinetics workflow](../../evaluation/kinetics_recovery/README.md)
defines the frozen inputs and rendering contract.
