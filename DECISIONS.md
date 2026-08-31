# Release decisions still requiring review

These items are deliberately unresolved in Phase 1 because each changes the public release contract.

## Repository and licensing

- **Public repository name and URL:** not chosen.
- **Code license:** not chosen. A license must be compatible with redistributed upstream method code and assets; method-specific licenses need review before copying code.
- **Citation metadata:** create `CITATION.cff` only after title, author order, repository URL, and archival DOI are stable.

## Data and large artifacts

- **Archival host:** decide whether processed benchmark inputs and selected method outputs will be deposited in Zenodo, Figshare, an institutional repository, or referenced from original sources.
- **Deposit granularity:** decide whether to publish only canonical processed inputs, also canonical final method outputs, or both.
- **Checksums and accession IDs:** add only after the deposit is frozen.

## Reproduction boundary

- **Method code:** decide case by case whether to redistribute a thin runner, vendor a pinned upstream snapshot, or provide an installation wrapper and upstream commit.
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
