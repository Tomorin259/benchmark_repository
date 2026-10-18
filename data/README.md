# Get the data

Choose inputs for method execution, or saved outputs for metric/figure
reproduction. Paths must match the selected module's input layout.

| Needed material | Access | Next step |
| --- | --- | --- |
| Canonical expression and simulation inputs | [Hugging Face inputs](https://huggingface.co/datasets/Zelin259/rna-velocity-benchmark-inputs), distribution `v1.0.0` | Select direct URLs in [input_file_manifest.tsv](input_file_manifest.tsv), save at their relative paths and check size/SHA-256 |
| Battich kinetic references, Tables S1/S3/S4 | Original study | Follow [UPSTREAM_REFERENCES.md](UPSTREAM_REFERENCES.md) and check the obtained CSV identities |
| Saved method outputs and evaluation assets | [Hugging Face results](https://huggingface.co/datasets/Zelin259/rna-velocity-benchmark-results) | Consult the actual availability index and [extraction guide](RESULT_DOWNLOADS.md) |
| Compact scores and plot-data | Owning `evaluation/` or `supplement/` module in Git | Follow the [evaluation index](../evaluation/README.md) or [figure map](../docs/figure_map.md) |

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

## File identities and layouts

| File | Purpose |
| --- | --- |
| [datasets.tsv](datasets.tsv) | Dataset ID, role, dimensions and upstream source |
| [input_file_manifest.tsv](input_file_manifest.tsv) | Direct URLs, paths, sizes and hashes for the 25 mirrored files |
| [input_release_manifest.tsv](input_release_manifest.tsv) | All 28 logical input dependencies and processing roles, including three upstream-only references |
| [upstream_references.tsv](upstream_references.tsv) | Original table identity, expected benchmark CSV hashes and supported layouts |
| [release_map.tsv](release_map.tsv) | Logical artifacts and materialization paths |
| [artifacts.tsv](../artifacts.tsv) | Artifact identity, checksum and registered access |

The [reference guide](UPSTREAM_REFERENCES.md) describes the expected benchmark
CSV format and checksums. [RESULT_DOWNLOADS.md](RESULT_DOWNLOADS.md) explains
archive/member hashes, logical aliases and the current availability index.

## Release preparation

For maintainers with the original artifacts, check a selected package:

```bash
python data/prepare_release.py --artifact-root /path/to/velocity-artifacts \
  --package reference-moments --output-root outputs/reference-release-check
```

Add `--build` to create a lossless TAR.GZ after source checks pass.
Publish the resulting package and register its download URL as a separate release step.
