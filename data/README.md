# Get the data

Choose inputs for method execution, or saved outputs for metric/figure
reproduction. Paths must match the selected module's input layout.

| Needed material | Access | Next step |
| --- | --- | --- |
| Canonical expression and simulation inputs | [Hugging Face inputs](https://huggingface.co/datasets/Zelin259/rna-velocity-benchmark/tree/v1.0.0/v1.0.0/inputs), distribution `v1.0.0` | Select direct URLs in [input_file_manifest.tsv](input_file_manifest.tsv), save at their relative paths and check size/SHA-256 |
| Battich kinetic references, Tables S1/S3/S4 | Original study | Follow [UPSTREAM_REFERENCES.md](UPSTREAM_REFERENCES.md) and check the obtained CSV identities |
| Saved method outputs and evaluation assets | [Hugging Face results v1.0.0](https://huggingface.co/datasets/Zelin259/rna-velocity-benchmark/tree/v1.0.0) | Select [download packages](result_package_manifest.tsv), map their [original members](result_member_manifest.tsv), and follow the [extraction guide](RESULT_DOWNLOADS.md) |
| Compact scores and plot-data | Owning `evaluation/` or `supplement/` module in Git | Follow the [evaluation index](../evaluation/README.md) or [figure map](../docs/figure_map.md) |
| Figure/table-indexed Source Data | [Source Data v1.0.0](SOURCE_DATA.md) | Browse individual tables and saved plot inputs; file sizes and hashes are in [source_data_manifest.tsv](source_data_manifest.tsv); a complete ZIP is optional |

## Canonical inputs

The input release has 25 original-format scientific files: 22 H5AD, one loom
and two gene-usage CSVs, covering 17 benchmark datasets. These are individual
downloads, with no TAR extraction. Dataset-specific source terms apply;
the release does not assign one new license to all third-party data.

The scEU labelled sources retain `sl`, `su`, `ul` and `uu` layers:
RPE1 has 5,422 cells × 11,848 genes, and organoid has 3,831 × 9,157.
They are separate from the benchmark main inputs and from Dynamo method outputs.
RPE1 alignment uses `Gene_Id`; the organoid main input lies within the labelled
source's cell/gene scope.

The [reference-gene usage tables](kinetics_reference_gene_usage/README.md)
describe the expanded 250-gene RPE1 and 270-gene organoid evaluation, selected
Figure 6 display flags, and valid/undefined/unavailable coverage. Original rate
trajectories are obtained through the upstream reference-table route.

## Local data roots

Choose the outer directories on your machine. The evaluators use these
different roots, according to their existing input maps:

| Command / parameter | Directory to supply | Example file below that directory |
| --- | --- | --- |
| Real-data method `--data-root` or top-level `DATA_ROOT` | Scientific input directory | `Pancreas/Pancreas.h5ad` |
| Simulation method `--data-root` | Velocity directory containing the registered simulation paths | `indicators/gene/simulation/branching_analyze/dataset/s1_branch.h5ad` |
| Cell and simulation-Cell `--artifact-root` | Parent of a `velocity/` directory | `velocity/results/...h5ad` |
| Gene, ensemble, veloVI coordinate diagnostic and gene-case `--artifact-root` | The velocity directory itself | `results/...h5ad` or `indicators/gene/...h5ad` |
| Cell `--reference-moments-root` | Directory of downloaded reference NPZs | `<evaluation_dataset>.npz` |
| Kinetics `--artifact-root` | Directory containing the reference and method aliases | `references/kinetics/rpe1_table_s1.csv`, `method-results/kinetics/...h5ad` |

For example, a Cell artifact root can contain `velocity/data/`,
`velocity/results/` and `velocity/indicators/`. Gene and ensemble receive
that root's `velocity/` subdirectory; real-method input roots can point to its
`velocity/data/` subdirectory. Preserve each input-file manifest path below
the velocity directory. Simulation paths begin with `indicators/`, and their
method commands receive the velocity directory itself.
The checkout's own `data/` directory holds the public registries and guides.

`velocity:<relative-path>` and `artifact://velocity/<relative-path>` identify
the same logical velocity-root file. A `repository:` locator is relative to
the checkout. Kinetics aliases are recorded alongside original paths in its
[run manifest](../evaluation/kinetics_recovery/results/central_results/manifest.json).
Map a downloaded member to the alias used by the selected evaluator, retaining
its bytes and checksum. Archive member names are transport paths; the
[download guide](RESULT_DOWNLOADS.md) explains their mapping.

## File identities and layouts

Dataset IDs follow `dataset_id` in [datasets.tsv](datasets.tsv); figures and
tables use the corresponding `paper_name`.
`release_status=mirrored_input` identifies inputs available in the public
manifest; the source/accession fields retain their original study attribution.

| File | Purpose |
| --- | --- |
| [datasets.tsv](datasets.tsv) | Dataset ID, role, dimensions and upstream source |
| [input_file_manifest.tsv](input_file_manifest.tsv) | Direct URLs, paths, sizes and hashes for the 25 mirrored files |
| [input_release_manifest.tsv](input_release_manifest.tsv) | All 28 logical input dependencies and processing roles, including three upstream-only references |
| [result_package_manifest.tsv](result_package_manifest.tsv) | 339 version-pinned result-package URLs, sizes and archive SHA-256 values |
| [result_member_manifest.tsv](result_member_manifest.tsv) | 388 original scientific files, their packages, member sizes/hashes, source locators and logical aliases |
| [upstream_references.tsv](upstream_references.tsv) | Original table identity, expected benchmark CSV hashes and supported layouts |
| [release_map.tsv](release_map.tsv) | Logical artifacts and materialization paths |
| [artifacts.tsv](../artifacts.tsv) | Artifact identity, checksum and registered access |

The [reference guide](UPSTREAM_REFERENCES.md) describes the expected benchmark
CSV format and checksums. [RESULT_DOWNLOADS.md](RESULT_DOWNLOADS.md) explains
archive/member hashes, logical aliases and the current availability index.
