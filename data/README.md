# Dataset metadata

`datasets.tsv` records the paper dataset identity, role, dimensions, stable
upstream source, and release status. It contains no server working paths.

Large canonical inputs are stored in the
[public Hugging Face input release](https://huggingface.co/datasets/Zelin259/rna-velocity-benchmark-inputs-public), distribution **v1.0.1**.
It covers 17 benchmark datasets and mirrors **25 scientific files** in 18
archives. Three original Battich reference tables are obtained from the
original study and are not mirrored; follow [UPSTREAM_REFERENCES.md](UPSTREAM_REFERENCES.md).
Full method-result H5AD downloads remain a separate, unpublished release.

## Input files, packages, and upstream references

- `input_release_manifest.tsv` inventories all 28 logical dependencies with
  roles, dimensions/layers and scientific-file SHA-256 values. Its storage
  fields distinguish 25 `huggingface_archive` records from three `upstream_only`
  records. An upstream-only row has no mirrored download URL.
- `input_package_manifest.tsv` maps the 25 bundled files to their 18 archives,
  public URLs, archive sizes/SHA-256, and member paths/sizes/SHA-256. A URL
  downloads the entire archive; it is not a direct H5AD member URL. Extract
  selected archives under your artifact root and verify the appropriate hash.
- `upstream_references.tsv` identifies Tables S1/S3/S4 by DOI and table number,
  records expected benchmark CSV identities and both supported input layouts.
  `check_references.py` checks locally supplied files without downloading or
  rewriting them. Current publisher downloads and any format conversion have
  not been verified as byte-identical to the recorded benchmark CSVs.
- `release_map.tsv` continues to map external logical artifacts to their
  materialization paths; it is not a list of archive download URLs.

The 23 expression files and two project-prepared gene-usage CSVs preserve
their previous scientific SHA-256 values. The 17 expression archives retain
their frozen bytes. `reference-gene-usage.tar.gz` contains only the two usage
tables, not the original Battich Tables S1/S3/S4. The old private reference
archive is not part of this public release.

Source attribution and dataset-specific terms are documented in the HF
release. The project does not relabel all third-party data under one new
license or assert that absence of a verified mirror license prohibits reuse
of numerical facts. Original studies retain attribution.

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
publication in the selected artifact store. Input public storage is linked above; other external artifact packages remain pending.
