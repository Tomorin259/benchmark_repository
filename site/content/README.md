# Website content sources

`website_content.json` contains reader-facing context, not evaluation values
or a new method/dataset registry. IDs are joined to `methods/configurations.tsv`
and `data/datasets.tsv`; missing or extra IDs stop export. Original scientific
scores, images and curve series retain their owning-module sources.

The content provenance records the SHA-256 of `merged.tex` read directly from
the online editor on 2026-10-13. Its source sections are recorded in JSON. That
read is an observed manuscript state, not a claim that subsequent paper edits
have been synchronized automatically. Figure numbering, captions and composites remain separate from browsing
existing dataset/method/diagnostic images.

- Dataset modality, context and trajectory terms come from the manuscript's
  dataset table. Dimensions/accessions stay owned by `data/datasets.tsv`.
  The scNT neuron benchmark subset (2,000 genes) is distinguished from the
  manuscript source-object context (44,021 genes).
- Model parameter assumptions and native `du/dt` descriptions are transcribed
  from the retained text table in `figures/assembly/fig7/fig7.pptx`.
  A dash is an absent entry in that table, not proof of biological incapacity.
  Gene evaluation coverage instead follows the saved spliced/unspliced layer
  contract in `evaluation/gene/config.py`. Repository post-processing can add
  a component not exported by a native model workflow.
- Cell-specific parameter capabilities follow the online kinetic-recovery
  text. Scored parameters are derived from the frozen result rows. These two
  properties are intentionally separate; beta is not a scored recovery metric.
  Labelled Dynamo is separate from the 20 RNA-seq configurations.
- Simulation dimensions are verified against the dataset registry and scenario
  config. Counts describe the generated object; evaluated special genes are a
  smaller saved selection. New generator invocations remain unseeded independent
  realizations.
- Ensemble values and ranks are copied from `ensemble/results/`. Dataset IDs
  use the Cell config map (with its documented dentate alias). Macro averages
  are the frozen metadata values; displayed differences are arithmetic deltas.
- Metric definitions, denominators and missing-value treatment follow the
  manuscript Metrics section and owning module guides. CosSim counts are not
  borrowed from CBDir. Organoid branch and system weighting remains unchanged.

Content changes require a source read and review; they must not update scientific
CSV files to make a statement true. Detailed audit downloads copy original bytes.
The exporter records the historical visual bundle's materialization commit and
its validation against the current build commit separately; a different commit
identity alone does not invalidate hash-identical display assets.

`diagnostic_assets.json` pins the existing diagnostic images selected for website
browsing. Dataset identities follow the owning veloAE renderer and projection
manifest; the PC1 summary uses its hash-verified SVG. TIvelo sources follow the
approved assembly and frozen paired-score data. The build validates these image
bytes before copying them, without regenerating scientific figures or changing
an owning manifest. The registration commit is not a new generation timestamp.
