# Dataset metadata

`datasets.tsv` records the paper dataset identity, role, dimensions, stable
upstream source, and release status. It contains no server working paths.

Large canonical inputs are not stored in Git. The root
[`artifacts.tsv`](../artifacts.tsv) records their logical identities, expected
checksums, and intended external-location fields. Stable public URIs have not
yet been published, so the registry is not currently a download index. A
method script may select a method-specific artifact, such as a `_counts.h5ad`
input, when its official workflow requires raw counts.

## External release preparation

`release_map.tsv` covers every external artifact ID exactly once. It adds only
package/member layout and logical source locators; hashes and publication state
remain owned by `artifacts.tsv`. Locators recovered from frozen Cell, kinetics,
veloVI and simulation manifests are not treated as current-byte verification.
Unfrozen input identity and unrecovered locators remain explicit.

## Input release inventory

`input_release_manifest.tsv` is the current release-preparation inventory for
the 28 source files selected for the 17 benchmark datasets: 14 benchmark main
inputs, four counts inputs, two scEU labelled sources, three fixed simulation
realizations, three Battich reference tables, and two expanded reference-gene
usage tables. It records dataset/role/package identity, member path, bytes,
SHA-256, H5AD dimensions/layers where applicable, and current publication
status.

The server-side file bytes and SHA-256 values were checked on 2026-10-12, but
the manifest intentionally has no private source paths or public URLs. The
`server_stat_sha256_verified_pending_publication` status means that the source
was found and hashed on the preparation host; it is not a Zenodo or Google
Drive publication claim. `public_storage=zenodo_core_pending` remains empty
until a real upload is completed and the download is verified.

The two scEU labelled sources are the native saved objects, not method-result
files: RPE1 retains `sl`, `su`, `ul`, and `uu` layers (5,422 cells × 11,848
genes), while organoid retains the same four labelled layers (3,831 cells ×
9,157 genes). The benchmark main inputs remain separate records with their
own roles and hashes. The RPE1 main input is aligned to the labelled source by
`Gene_Id`; the organoid main input is within the labelled source's cell/gene
scope. No source file is replaced by a Dynamo output.

`kinetics_reference_gene_usage/` contains the 250-gene RPE1 and 270-gene
organoid tables, their manifest, and an explanation of valid/undefined/
unavailable status. These tables describe the expanded reference evaluation;
they do not duplicate the Table S1/S3/S4 trajectories or claim to describe
model-training gene usage.

Check one package without creating a large archive:

```bash
python data/prepare_release.py --artifact-root /path/to/velocity-artifacts \
  --package reference-moments --output-root outputs/reference-release-check
```

Add `--build` to write a `.tar.gz` only after every selected source matches its
registered SHA-256. The archive preserves the logical velocity-root paths used
by the existing launchers/evaluators. Missing sources or changed bytes block
archive construction; no new canonical input or score is inferred. The script
never uploads or fills download URIs. Record real URIs only after an actual
publication in the selected artifact store. Public storage selection is pending.
