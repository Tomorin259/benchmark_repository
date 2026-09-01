# Release decisions still requiring review

These items remain deliberately unresolved after the inventory and execution-registry phases because each changes the public release contract.

## Repository and licensing

- **Public repository name and URL:** not chosen.
- **Code license:** not chosen. A license must be compatible with redistributed upstream method code and assets; method-specific licenses need review before copying code.
- **Citation metadata:** create `CITATION.cff` only after title, author order, repository URL, and archival DOI are stable.

## Data and large artifacts

- **Canonical processed-input archive:** unresolved; freeze checksum-addressed inputs.
- **Derived-input preprocessing reconstruction:** unresolved for transformed or gene-subset H5ADs.
- **Simulation regeneration:** unresolved until seeds, output paths, and revisions are recorded.
- **Metric row-level input manifest:** unresolved; current entrypoint has empty `ADATA_PATHS`.
- **DG_lamanno variant handling:** unresolved; preserve both H5AD variants pending review.
- **Upstream source redistribution/licensing:** unresolved for processed inputs and kinetic reference tables.

- **Archival host:** decide whether processed benchmark inputs and selected method outputs will be deposited in Zenodo, Figshare, an institutional repository, or referenced from original sources.
- **Deposit granularity:** decide whether to publish only canonical processed inputs, also canonical final method outputs, or both.
- **Checksums and accession IDs:** add only after the deposit is frozen.

## Reproduction boundary

- **Method code:** decide case by case whether to redistribute a thin runner, vendor a pinned upstream snapshot, or provide an installation wrapper and upstream commit.
- **Runner promotion:** curate one method-specific script at a time with an explicit dataset selector or mapping, output contract, and reviewed environment specification. A shared universal runner is not required, and a documented manual copy or rename of final outputs is acceptable.
- **TFvelo release curation:** the source tree is located at `project:TFvelo-main`, and the 14 canonical result files match the original result tree by SHA-256. The remaining decisions are the thin-runner interface, upstream revision, and redistribution terms.
- **veloVI release contract:** the cell-level `without_vu` family is confirmed as the default spliced `get_velocity` output, while the gene-level family additionally stores the explicit unspliced output as `velocity_u`. The remaining work is to curate the two runner roles and their clean environment specification; there is no unresolved output-mode conflict.
- **TSvelo release curation:** the current render-manifest paths are the canonical 14 outputs; `*_ori`, `*_1`, and `results_old` are historical output variants. The remaining decisions are the preprocessing/input mapping, method revision, and bundled-asset license terms.
- **Historical environments:** resolve ambiguous or conflicting assignments and pin upstream revisions before writing clean environment specifications; raw Conda exports are not accepted as release files.
- **Execution tiers:** decide which workflows are required as smoke tests and which remain full HPC reproductions.
- **Containers:** decide whether per-method Conda environments are sufficient or whether selected workflows also need container images.

## Figures

- **Figure 1:** no editable source or complete generator was found in the working project. Decide whether to release only the final PDF or recover the original editable asset.
- **Figure 4 and Figures S1–S8:** canonical result collections are known, but final panel assembly is partly manual. Panel-level source packages must be frozen during the figure-reproduction phase.
- **Figure 6:** direct summary CSV files are known, but no complete script that rebuilds the entire final figure has been found. Decide whether to implement a clean reproduction script or document manual assembly.
- **Figures S9–S12:** the current paper builder is in the paper snapshot and must be migrated together with its small inputs during the ensemble phase.

## Publication text

- **Data Availability:** replace the repository placeholder only after repository and archive URLs are persistent.
- **Code Availability:** state clearly which stages are fully reproducible, which require external method repositories, and which final figures include documented manual assembly.
