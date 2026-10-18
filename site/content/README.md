# Website content sources

`website_content.json` contains reader-facing context. IDs are joined to `methods/configurations.tsv`
and `data/datasets.tsv`; missing or extra IDs stop export. Original scientific
scores, images and curve series retain their owning-module sources.

The content provenance records the SHA-256 of `merged.tex` read directly from
the online editor on 2026-10-13, with source sections recorded in JSON. Refresh
the relevant website text when the manuscript changes. Figure numbering,
captions and composite bindings follow the paper figure map.

- Dataset modality, context and trajectory terms come from the manuscript's
  dataset table. Dimensions/accessions stay owned by `data/datasets.tsv`.
  The scNT neuron benchmark subset (2,000 genes) is distinguished from the
  manuscript source-object context (44,021 genes).
- Model parameter assumptions and native `du/dt` descriptions are transcribed
  from the retained text table in `figures/assembly/fig7/fig7.pptx`.
  A dash denotes an unspecified entry in that source table.
  Gene evaluation coverage instead follows the saved spliced/unspliced layer
  contract in `evaluation/gene/config.py`. Repository post-processing can add
  a component not exported by a native model workflow.
- Cell-specific parameter capabilities follow the online kinetic-recovery
  text. The frozen result rows identify the scored alpha and gamma profiles.
  Labelled Dynamo is separate from the 20 RNA-seq configurations.
- Simulation dimensions are verified against the dataset registry and scenario
  config. Counts describe the generated object; evaluated special genes are a
  smaller saved selection. The generator creates additional simulated datasets.
- Ensemble values and ranks are copied from `ensemble/results/`. Dataset IDs
  use the Cell config map (with its documented dentate alias). Macro averages
  are the frozen metadata values; displayed differences are arithmetic deltas.
- Metric definitions, denominators and missing-value treatment follow the
  manuscript Metrics section and owning module guides. CosSim and CBDir use
  their respective valid-cell definitions. Organoid summaries retain the
  recorded branch and system weighting.

Review content changes against their scientific sources. Audit downloads copy
the source files. The exporter records the visual bundle's materialization
commit, the current validation commit and the checked asset hashes.

`diagnostic_assets.json` pins the existing diagnostic images selected for website
browsing. Dataset identities follow the owning veloAE renderer and projection
manifest; the PC1 summary uses its hash-verified SVG. TIvelo sources follow the
approved assembly and frozen paired-score data. The build validates these image
bytes before copying them. The owning manifests retain the original generation
identity of each asset.
