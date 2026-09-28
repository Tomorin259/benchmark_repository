# Paper figures

Paper-level PDFs and assembly assets are tracked here. The assembly manifest
maps each materialized panel to a repository source and SHA256 and records
unresolved independent sources without inventing replacements.

The online-matched PDFs are `figThree.pdf` (selected gene/cell transitions),
`figFour.pdf` (simulation), and `figFive.pdf` (real-data gene kinetics).
`figSix.pdf` is the validated user-exported scEU-seq kinetic-parameter recovery
composite, including the current B-panel palette and F-panel four-decimal means.
Figure 6 has 16 source panels, including the Dynamo(label data) organoid E
stream; `fig6_composite_manifest.json` records their hashes, the final PDF
checksum, and the verification scope. `figSeven.pdf` retains the repository-baseline integrated summary
values. The historical kinetics PDF bearing the old Figure 5 filename is
retained in `archive/kinetics/figFive.pdf`.

Do not add intermediate panels or alternate result copies to the active
figure directory.
