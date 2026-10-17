# Paper figures

`suppFigFifteen.pdf` is the retained author-exported MultiVelo diagnostic
composite. `suppfig15_composite_manifest.json` maps its A/B materials,
editable PPTX, compact numerical sources and Supplementary Table S9. Source
panels are owned by `supplement/multivelo/`; the assembly entry point under
`assembly/suppfig15/` preserves source scale and leaves final PPTX layout and
PDF export to the author.

Paper-level PDFs and assembly assets are tracked here. The assembly manifest
maps each referenced panel to a repository source and SHA256 and records
unresolved independent sources without inventing replacements.

`figOne.pdf`, `suppFigNine.pdf`, and `suppFigEleven.pdf` are retained original
PDFs whose git-blob hashes match the current online project. Their composite
manifests distinguish retained final exports from independent reconstruction.

The active PDFs are `figThree.pdf` (selected gene/cell transitions),
`figFour.pdf` (simulation), and `figFive.pdf` (real-data gene kinetics).
Figure 4 is the approved user-exported, cropped composite containing the paired
Cell/Gene D panel. `fig4_composite_manifest.json` records its checksum, source
panel, and manual-assembly boundary; online byte identity was not rechecked.
`figSix.pdf` is the validated user-exported scEU-seq kinetic-parameter recovery
composite, including the current B-panel palette and F-panel four-decimal means.
It uses `Dynamo (scEU)` display labels, two-decimal C-panel `Mean (n)` summary
cells sharing the gene heatmap's color scale and row boundaries, and D-panel
category headings with Macro-F1 shown below the Venn diagrams (0.324 for
cellDancer and 0.340 for DeepVelo). Cosine intervals are documented in the
kinetics module rather than repeated in D. C retains its source aspect ratio;
D--F use the approved compact layout. Figure 6 has 17 source panels, including
the Dynamo(label data) and TSvelo organoid E streams. A and E follow the same
five-method display order. `fig6_composite_manifest.json` records all source
hashes, the online-matched final PDF checksum, and the verification scope.
`figTwo.pdf`, `figSeven.pdf` and `suppFigTen.pdf` are the verified cropped author
exports containing the corrected Forebrain Cell CBDir and shared-frame Ensemble
BATC results. Their composite manifests record source and export hashes and the
crop boxes. The editable sources are `assembly/fig2/fig2.pptx`,
`assembly/fig7/fig7.pptx` and `assembly/suppfig10/suppfig10.pptx`.
The historical kinetics PDF bearing the old Figure 5 filename is retained in `archive/kinetics/figFive.pdf`.

Do not add intermediate panels or alternate result copies to the active
figure directory.

`suppFigTwelve.pdf` is the approved user-exported, cropped simulation Cell
composite. `suppfig12_composite_manifest.json` records its checksum and
manual-assembly boundary; numerical results and the Supplementary Table 7
exporter belong to `evaluation/cell/simulation/`.
