# Native result downloads

The [Hugging Face output repository](https://huggingface.co/datasets/Zelin259/rna-velocity-benchmark-results)
publishes files in verified batches. Its
[availability index](https://huggingface.co/datasets/Zelin259/rna-velocity-benchmark-results/blob/main/available_files.tsv)
lists committed and verified packages; the
[upload status](https://huggingface.co/datasets/Zelin259/rna-velocity-benchmark-results/blob/main/UPLOAD_STATUS.json)
records progress toward the complete collection.

The frozen selection contains 384 original scientific files: 370 H5AD and
14 NPZ, totaling 252,148,325,990 source bytes. Distribution uses 324 lossless
TAR.GZ archives and 14 unchanged NPZ files, totaling 149,558,222,960 download
bytes. These are method outputs and evaluation assets; original-format
benchmark inputs remain in the separate
[input release](https://huggingface.co/datasets/Zelin259/rna-velocity-benchmark-inputs).

The output repository provides:

- `available_files.tsv`: actual committed objects, sizes, SHA-256 values and
  immutable download URLs. Download only the packages you need.
- `v1.0.0/metadata/package_manifest.tsv`: every selected original file, its
  package, member path, original size/SHA, logical aliases and source locator.
- `v1.0.0/metadata/SHA256SUMS.txt`: compressed-package and metadata hashes.

Check the downloaded archive with `sha256sum` (Linux) or `shasum -a 256`
(macOS), then extract each archive into a separate directory:

```bash
mkdir -p unpacked/capo_pre_C09
tar -xzf C09.tar.gz -C unpacked/capo_pre_C09
```

For this package the original file is
`unpacked/capo_pre_C09/cell/cell2fate_capo_pre.h5ad`. Check its size and
SHA-256 against the member columns of the package manifest. **Archive hashes
and original-file hashes identify different byte streams**; native artifact
hashes in `artifacts.tsv` are not replaced by compressed-package hashes.
Large members use standard POSIX PAX extensions. NPZ files need no TAR
extraction. Allow space for both the archive and extracted data.

Map extracted members to the evaluator's expected logical source paths using
`source_locator`, registered artifact IDs and the owning module README.
Arrange the extracted files under the evaluator's input root. Files shared by
several evaluations are stored once, with their logical aliases retained.

The selection is frozen at code commit
`0e38d4fef12f71d8e879088de1b15e5510cc40a5`. It covers Cell outputs,
additional Gene/Kinetics outputs, fixed simulation outputs, ensemble outputs,
TIvelo controls and evaluation reference moments. Later diagnostic previews,
checkpoints, logs, environments and failed runs are outside this frozen set.
Packaging preserves the original scientific-file bytes. Compact scores and
figures remain in their owning GitHub modules.

Retain original dataset and method citations and dataset-specific terms; see
the output repository's SOURCE_TERMS.md and CITATION.cff. Three original
Battich kinetic-reference tables are obtained from the original study through
the [reference access guide](UPSTREAM_REFERENCES.md).
