# Paper figure and table map

Find a final figure first, then follow its source module for numerical data and
reproduction commands. Final composites and experiment-level source panels
have different reconstruction requirements.

## Main figures

| Figure | Final PDF | Numerical sources and commands | Assembly |
| --- | --- | --- | --- |
| 1 — benchmark framework | [figOne.pdf](../figures/figOne.pdf) | Study-design schematic | [Retained author source](../figures/assembly/fig1/authoring/source.pptx); independent sources for the corrected final panels are incomplete |
| 2 — Cell benchmark | [figTwo.pdf](../figures/figTwo.pdf) | [Cell tables and plots](../evaluation/cell/README.md) | [Editable assembly](../figures/assembly/fig2/fig2.pptx); author PDF export |
| 3 — selected gene–transition cases | [figThree.pdf](../figures/figThree.pdf) | [Selected-edge panels](../evaluation/gene/figures/panels/selected_edges/) and [Gene guide](../evaluation/gene/README.md) | Retained composite; three panel sources are registered in the assembly map |
| 4 — simulated velocity evaluation | [figFour.pdf](../figures/figFour.pdf) | [Gene panels](../evaluation/gene/README.md); final D uses [paired Cell/Gene results](../evaluation/cell/simulation/README.md) | [Author source](../figures/assembly/fig4/authoring/source.pptx); author PDF export |
| 5 — real-data Gene evaluation | [figFive.pdf](../figures/figFive.pdf) | [Gene results and panels](../evaluation/gene/README.md) | [Author source](../figures/assembly/fig5/authoring/source.pptx); D overview sources follow the assembly manifest |
| 6 — scEU-seq parameter recovery | [figSix.pdf](../figures/figSix.pdf) | [Kinetics guide](../evaluation/kinetics_recovery/README.md); [exact plotted rows](../evaluation/kinetics_recovery/figures/fig6/plot_data/) | Quantitative panels are scripted; labelled Dynamo A/E streams retain manual assets |
| 7 — integrated comparison | [figSeven.pdf](../figures/figSeven.pdf) | [Summary tables and refresh commands](../evaluation/summary/README.md) | [Editable charts](../figures/assembly/fig7/fig7.pptx); manual WPS export and crop |

## Supplementary figures

For S1–S8 and S14, the final composite is retained in the manuscript project;
the links below lead to their repository source panels and numerical data.

| Figure | Subject | Sources and reproduction |
| --- | --- | --- |
| S1 | Dentategyrus_lamanno branching genes | [Gene panels](../evaluation/gene/figures/panels/real/DentateGyrus_lamanno/), [renderer](../evaluation/gene/README.md) |
| S2 | Mouse retina branching genes | [Gene panels](../evaluation/gene/figures/panels/real/retina/), [renderer](../evaluation/gene/README.md) |
| S3 | Gastrulation MURK genes, part I | [Gene panels](../evaluation/gene/figures/panels/real/gastrulation/), [renderer](../evaluation/gene/README.md) |
| S4 | Gastrulation MURK genes, part II | [Gene panels](../evaluation/gene/figures/panels/real/gastrulation/), [renderer](../evaluation/gene/README.md) |
| S5 | capofibro MURK genes | [Gene panels](../evaluation/gene/figures/panels/real/capo_pre/), [renderer](../evaluation/gene/README.md) |
| S6 | Human bone marrow MURK genes | [Gene panels](../evaluation/gene/figures/panels/real/bonemarrow/), [renderer](../evaluation/gene/README.md) |
| S7 | mousebrain MURK genes | [Gene panels](../evaluation/gene/figures/panels/real/mousebrain/), [renderer](../evaluation/gene/README.md) |
| S8 | Pancreas MURK genes | [Gene panels](../evaluation/gene/figures/panels/real/pancreas/), [renderer](../evaluation/gene/README.md) |
| S9 | Matched Cell/Gene comparison | [Retained PDF](../figures/suppFigNine.pdf); [matched Cell comparison](../evaluation/cell/README.md#gene-output-matched-cell-comparison) and [Gene scores](../evaluation/gene/README.md); independent final-composite reconstruction remains incomplete |
| S10 | GeoMedian-4 Cell evaluation | [Retained PDF](../figures/suppFigTen.pdf); [ensemble renderer](../ensemble/README.md), [editable assembly](../figures/assembly/suppfig10/suppfig10.pptx) |
| S11 | veloAE representation diagnostics | [Retained PDF](../figures/suppFigEleven.pdf); [PCA inputs and renderers](../supplement/veloae/README.md); final whole-page reconstruction remains incomplete |
| S12 | Cell evaluation on three simulations | [Retained PDF](../figures/suppFigTwelve.pdf); [scores, plots and table](../evaluation/cell/simulation/README.md) |
| S13 | TIvelo root and ordering controls | [Retained PDF](../figures/assembly/suppfig13/suppFigThirteen.pdf); [diagnostics](../supplement/tivelo_gene/README.md), [editable assembly](../figures/assembly/suppfig13/README.md) |
| S14 | Expanded reference-gene parameter recovery | [Expanded scores](../evaluation/kinetics_recovery/results/expanded_reference/README.md), [renderers](../evaluation/kinetics_recovery/plotting/README.md) |
| S15 | RNA-only MultiVelo phase constraints | [Retained PDF](../figures/suppFigFifteen.pdf); [controls and table](../supplement/multivelo/README.md), [assembly](../figures/assembly/suppfig15/) |

## Supplementary tables

| Table | Subject | Repository source |
| --- | --- | --- |
| S1 | Method versions and execution settings | [Method registry and launchers](../methods/README.md), [environment limits](environments.md) |
| S2 | Dataset descriptions | [Dataset registry](../data/datasets.tsv), [simulation configuration](../evaluation/gene/simulation/config.json) |
| S3 | Cell scores, continued across parts | [Primary Cell table](../evaluation/cell/results/cell_metrics.csv) |
| S4 | Cross-dataset Cell summaries | [Summary guide](../evaluation/summary/README.md) |
| S5 | Signed-transition ensemble | [Scores and independent ranking](../ensemble/README.md) |
| S6 | veloVI coordinate sensitivity | [CSV, transformation and table export](../supplement/velovi/README.md) |
| S7 | Cell scores on simulations | [Simulation Cell table exporter](../evaluation/cell/simulation/README.md) |
| S8 | Gene sets, coverage and scores | [Selected genes and screening evidence](../evaluation/gene/selection/README.md), [Gene scores](../evaluation/gene/results/README.md) |
| S9 | scVelo/MultiVelo phase controls | [Full-precision summary and LaTeX table](../supplement/multivelo/README.md) |

## Source identity and manual assembly

[assembly_manifest.tsv](../figures/assembly/assembly_manifest.tsv) gives exact
panel paths, checksums and missing-source flags. Figure-specific composite
manifests in [figures/](../figures/) record final export identity and its
verification scope; [the Gene panel manifest](../evaluation/gene/figures/panels/assembly_manifest.tsv)
records Gene panel/score bindings.

Use [material export](../figures/assembly/README.md) for fresh source copies.
After updating a source panel, assemble and export the corresponding composite.
Figure-specific source availability is recorded in the table and manifests above.
Historical kinetics assets named Figure 5 belong to the
former numbering; current parameter recovery is Figure 6.
