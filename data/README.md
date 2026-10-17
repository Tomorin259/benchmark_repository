# Dataset metadata

`datasets.tsv` records the paper dataset identity, role, dimensions, stable
upstream source, and release status. It contains no server working paths.

Large canonical inputs are stored as individual files in the
[public Hugging Face input release](https://huggingface.co/datasets/Zelin259/rna-velocity-benchmark-inputs), distribution **v1.0.2**.
It covers 17 benchmark datasets and provides **25 original-format scientific
files: 22 H5AD, one loom and two CSVs**. Select the files needed for your
analysis; no TAR extraction is required. Three original Battich reference
tables are obtained from the original study; follow
[UPSTREAM_REFERENCES.md](UPSTREAM_REFERENCES.md).

## Input files and upstream references

- `input_release_manifest.tsv` inventories all 28 logical dependencies with
  processing roles, dimensions/layers and scientific-file SHA-256 values.
  It distinguishes 25 `huggingface_file` records from three `upstream_only`
  records. An upstream-only row has no mirrored download URL.
- `input_file_manifest.tsv` maps the 25 individual files to direct public
  URLs, relative paths, sizes and SHA-256 values. Save a file at the recorded
  path under your artifact root and check its scientific-file hash.
- `upstream_references.tsv` records Tables S1/S3/S4 by DOI and table number,
  expected benchmark CSV identities and both supported input layouts.
  `check_references.py` checks locally supplied files without downloading or
  rewriting them. Current publisher downloads and any format conversion have
  not been verified as byte-identical to the recorded benchmark CSVs.
- `release_map.tsv` continues to map external logical artifacts to their
  materialization paths. It is distinct from the input download manifest.

The 23 expression files and two project-prepared gene-usage CSVs preserve
their frozen scientific bytes and SHA-256 values. Working inputs,
count-preserving supplements, native labelled sources and fixed simulations
keep their documented processing roles. The project does not relabel all
third-party data under one new license; source attribution and dataset-specific
terms are documented in the HF release.

Saved method outputs and evaluation assets are being uploaded to the separate
[Hugging Face output release](https://huggingface.co/datasets/Zelin259/rna-velocity-benchmark-results)
in verified batches. It is not yet a complete release. Use its actual
availability index and follow [RESULT_DOWNLOADS.md](RESULT_DOWNLOADS.md) for
TAR.GZ extraction and the separate archive/member SHA-256 checks.

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
