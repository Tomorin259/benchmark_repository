# Paper figures

Paper-level PDFs and assembly assets are tracked here. The assembly manifest
maps each materialized panel to a repository source and SHA256 and records
unresolved independent sources without inventing replacements.

The current online-matched PDFs are `figThree.pdf` (selected gene/cell
transitions), `figFour.pdf` (simulation), `figFive.pdf` (real-data gene
kinetics), and `figSix.pdf` (scEU-seq kinetic-parameter recovery). Figure 6
has 16 source panels, including the Dynamo(label data) organoid E stream;
`fig6_composite_manifest.json` records their hashes and the online PDF
checksum. `figSeven.pdf` retains the repository-baseline integrated summary
values. The historical kinetics PDF bearing the old Figure 5 filename is
retained in `archive/kinetics/figFive.pdf`.

Do not add intermediate panels or alternate result copies to the active
figure directory.
