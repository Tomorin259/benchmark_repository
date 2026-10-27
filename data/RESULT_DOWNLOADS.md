# Native result downloads

[Hugging Face v1.0.0](https://huggingface.co/datasets/Zelin259/rna-velocity-benchmark/tree/v1.0.0)
provides 339 download objects containing 388 scientific files: 374 H5ADs and
14 NPZs. The 325 TAR.GZ archives and 14 NPZ files total 150.61 GB.
[Benchmark inputs](https://huggingface.co/datasets/Zelin259/rna-velocity-benchmark/tree/v1.0.0/v1.0.0/inputs)
are available as individual files.

- [Package manifest](result_package_manifest.tsv): download URLs, sizes and archive SHA-256.
- [Member manifest](result_member_manifest.tsv): original file paths, sizes, hashes,
  evaluation roles, artifact IDs and source locators.
- [Source identities](source_provenance.tsv): native method fits and simulated observations.

The collection covers Cell outputs, additional Gene/Kinetics outputs, simulations,
ensembles, TIvelo/MultiVelo controls and reference moments. MultiVelo Cell files
include native fitted velocities, normalized spliced velocity, its graph and
projection. Compact scores and figures are in the owning repository modules.

## Download and extract

Select a package and verify its SHA-256 with `sha256sum` (Linux) or
`shasum -a 256` (macOS). Extract each archive into its own directory:

```bash
mkdir -p unpacked/capo_pre_C09
tar -xzf C09.tar.gz -C unpacked/capo_pre_C09
```

For this example, the member is
`unpacked/capo_pre_C09/cell/cell2fate_capo_pre.h5ad`. Verify it against the
member manifest. Archive hashes identify compressed bytes; member hashes
identify the original H5AD/NPZ bytes. NPZ files can be loaded directly with
`numpy.load(path, allow_pickle=False)`.

Allow space for both archive and extracted files. The largest archive is
18.81 GB and its H5AD is 32.13 GB. TAR.GZ uses standard POSIX PAX extensions.

## Arrange evaluation inputs

Use `source_locator` to place extracted members under the selected evaluator's
input root. For example,
`velocity:results/cell2fate/capo_pre/cell2fate_capo_pre.h5ad` maps to
`<velocity-root>/results/cell2fate/capo_pre/cell2fate_capo_pre.h5ad`.
The [local-root table](README.md#local-data-roots) defines Cell, Gene, ensemble
and kinetics layouts. Files shared across evaluations have multiple logical
aliases in the member manifest.

Hugging Face also provides `available_files.tsv` for all download objects,
`v1.0.0/metadata/SHA256SUMS.txt` and the package/member manifests.
Dataset and method attribution is recorded in `SOURCE_TERMS.md` and
`CITATION.cff`. Obtain the Battich kinetic-reference tables through the
[reference access guide](UPSTREAM_REFERENCES.md).
