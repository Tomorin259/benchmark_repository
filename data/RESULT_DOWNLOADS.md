# Native result downloads

The [Hugging Face results v1.0.1](https://huggingface.co/datasets/Zelin259/rna-velocity-benchmark/tree/v1.0.1)
contain 339 download objects, including four native MultiVelo candidate-gate
H5ADs in one 17.97 MB archive. The diagnostic archive and all four extracted
members were verified by anonymous download.

The index contains 388 original scientific files: 374 H5AD and
14 NPZ, totaling 253,295,202,347 source bytes. Distribution uses 325 lossless
TAR.GZ archives and 14 NPZ files, totaling 150,607,413,075 download bytes.
Original-format benchmark inputs remain in the separate
[input release](https://huggingface.co/datasets/Zelin259/rna-velocity-benchmark/tree/v1.0.1/v1.0.1/inputs).

Choose downloads and their original files through:

- [result_package_manifest.tsv](result_package_manifest.tsv): 339 package URLs
  pinned to `v1.0.1`, sizes and package SHA-256 values.
- [result_member_manifest.tsv](result_member_manifest.tsv): 388 original
  scientific files, their package IDs, archive member paths, original hashes,
  evaluation roles, registered aliases and source locators.

The output repository also provides:

- `available_files.tsv`: actual committed objects, sizes, SHA-256 values and
  immutable download URLs. Download only the packages you need.
- `v1.0.1/metadata/package_manifest.tsv`: every selected original file, its
  package, member path, original size/SHA, logical aliases and source locator.
- `v1.0.1/metadata/SHA256SUMS.txt`: compressed-package and metadata hashes.
- `v1.0.1/metadata/result_package_manifest.tsv` and
  `v1.0.1/metadata/result_member_manifest.tsv`: package and member download indexes.

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
The largest archive is about 18.81 GB and contains an approximately 32.13 GB
H5AD; allow room for both while verifying it.

Map extracted members to the evaluator's expected logical source paths using
`source_locator`, registered artifact IDs and the owning module README.
Arrange the extracted files under the evaluator's input root. Files shared by
several evaluations are stored once, with their logical aliases retained.

For example, the C09 member above has source locator
`velocity:results/cell2fate/capo_pre/cell2fate_capo_pre.h5ad`. Its evaluation
location is `<velocity-root>/results/cell2fate/capo_pre/cell2fate_capo_pre.h5ad`,
regardless of the archive's `cell/` directory. Cell receives the parent of this
velocity root. Gene and ensemble receive the velocity root directly for the
methods selected by their own input maps. The
[local-root table](README.md#local-data-roots) also covers method inputs,
reference moments and kinetics aliases. Copies or local links can expose
the required aliases while preserving the original-file identity.

Join a registered artifact ID to `registered_artifact_ids` in the member
manifest, then join `package_id` to the package manifest. The native artifact
hash identifies the original member; its URI stays empty when that original
file is transported inside a TAR.GZ rather than served individually.

The package and member manifests above define the published file selection.
It covers Cell outputs, additional Gene/Kinetics outputs, fixed simulation
outputs, ensemble outputs, TIvelo and MultiVelo controls and evaluation
reference moments.
Each package preserves the scientific-file bytes identified in its member manifest. Compact scores and
figures remain in their owning GitHub modules.

Retain original dataset and method citations and dataset-specific terms; see
the output repository's SOURCE_TERMS.md and CITATION.cff. Three original
Battich kinetic-reference tables are obtained from the original study through
the [reference access guide](UPSTREAM_REFERENCES.md).

MultiVelo Cell files contain native fitted velocities, the official normalized spliced velocity, and its matching graph and projection. Native fit and simulation source identities are listed in [source_provenance.tsv](source_provenance.tsv). Source locators define the evaluation aliases; their source-file hashes and package hashes remain separate.
