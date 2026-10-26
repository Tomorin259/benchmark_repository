# Source Data

[Browse Source Data v1.0.0](https://huggingface.co/datasets/Zelin259/rna-velocity-benchmark/tree/8c2d44d67c7113ad70707c415f13a58930079ec7/source_data/v1.0.0) for the main and supplementary figures
and tables. Files are available individually: 88 numerical/metadata tables,
243 saved NPZ plot inputs, and the Figure 1 framework schematic PDF.
Each figure/table directory describes its files and their corresponding panels.

| Material | Browse |
| --- | --- |
| Main figures | [Figures 1–7](https://huggingface.co/datasets/Zelin259/rna-velocity-benchmark/tree/8c2d44d67c7113ad70707c415f13a58930079ec7/source_data/v1.0.0/main_figures), including the Figure 1 schematic |
| Supplementary figures | [Figures S1–S15](https://huggingface.co/datasets/Zelin259/rna-velocity-benchmark/tree/8c2d44d67c7113ad70707c415f13a58930079ec7/source_data/v1.0.0/supplementary_figures) |
| Supplementary tables | [Tables S1–S9](https://huggingface.co/datasets/Zelin259/rna-velocity-benchmark/tree/8c2d44d67c7113ad70707c415f13a58930079ec7/source_data/v1.0.0/supplementary_tables) |
| Figure/panel-to-file index | [SOURCE_INDEX.tsv](https://huggingface.co/datasets/Zelin259/rna-velocity-benchmark/resolve/8c2d44d67c7113ad70707c415f13a58930079ec7/source_data/v1.0.0/SOURCE_INDEX.tsv?download=true) |
| Array field descriptions | [ARRAY_DICTIONARY.tsv](https://huggingface.co/datasets/Zelin259/rna-velocity-benchmark/resolve/8c2d44d67c7113ad70707c415f13a58930079ec7/source_data/v1.0.0/metadata/ARRAY_DICTIONARY.tsv?download=true) |

Download only the CSV, TSV or NPZ files you need; individual files require no
ZIP extraction. For all files, use the optional
[Source_Data_complete.zip](https://huggingface.co/datasets/Zelin259/rna-velocity-benchmark/resolve/8c2d44d67c7113ad70707c415f13a58930079ec7/source_data/v1.0.0/Source_Data_complete.zip?download=true) (174.94 MB).
It extracts into Source_Data/ with the same directory layout. The
[download manifest](source_data_manifest.tsv) records every individual file
and the optional archive, with immutable URLs, sizes and SHA-256 checksums.

You can also download the directory using the Hugging Face Python client:

```python
from huggingface_hub import snapshot_download

snapshot_download(
    repo_id="Zelin259/rna-velocity-benchmark",
    repo_type="dataset",
    revision="8c2d44d67c7113ad70707c415f13a58930079ec7",
    allow_patterns="source_data/v1.0.0/**",
    ignore_patterns="*.zip",
    local_dir="benchmark-results",
)
```

CSV/TSV files retain their saved precision. Load NPZ files with
`numpy.load(path, allow_pickle=False)` and consult the array field descriptions.
`metadata/` contains the array dictionary, phase-panel gene availability,
file checksums and a compact release record. Figure/table coverage and source
identities are recorded in `SOURCE_INDEX.tsv` and the figure/table directories.
Plotting code is maintained in the owning evaluation module. All five Figure 6A
methods use [the same scVelo stream entry](../evaluation/kinetics_recovery/plotting/README.md).

Native method outputs are available through the [result download guide](RESULT_DOWNLOADS.md),
[package index](result_package_manifest.tsv) and [member index](result_member_manifest.tsv).
