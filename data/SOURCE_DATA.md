# Source Data

[Browse Source Data v1.0.0](https://huggingface.co/datasets/Zelin259/rna-velocity-benchmark/tree/cd1adffe90b778ec1d11ed5985284f1a5d41095c/source_data/v1.0.0) for the main and supplementary figures
and tables. Files are available individually: 88 numerical/metadata tables,
243 saved NPZ plot inputs, plotting scripts and retained image materials.
Each figure/table directory describes its files and their corresponding panels.

| Material | Browse |
| --- | --- |
| Main figures | [Figures 2–7](https://huggingface.co/datasets/Zelin259/rna-velocity-benchmark/tree/cd1adffe90b778ec1d11ed5985284f1a5d41095c/source_data/v1.0.0/main_figures); Figure 1 is a schematic |
| Supplementary figures | [Figures S1–S15](https://huggingface.co/datasets/Zelin259/rna-velocity-benchmark/tree/cd1adffe90b778ec1d11ed5985284f1a5d41095c/source_data/v1.0.0/supplementary_figures) |
| Supplementary tables | [Tables S1–S9](https://huggingface.co/datasets/Zelin259/rna-velocity-benchmark/tree/cd1adffe90b778ec1d11ed5985284f1a5d41095c/source_data/v1.0.0/supplementary_tables) |
| Reproduction scripts and retained images | [Figure 5](https://huggingface.co/datasets/Zelin259/rna-velocity-benchmark/tree/cd1adffe90b778ec1d11ed5985284f1a5d41095c/source_data/v1.0.0/reproduction/fig5), [Figure 6A](https://huggingface.co/datasets/Zelin259/rna-velocity-benchmark/tree/cd1adffe90b778ec1d11ed5985284f1a5d41095c/source_data/v1.0.0/reproduction/fig6a) |
| Figure/panel-to-file index | [SOURCE_INDEX.tsv](https://huggingface.co/datasets/Zelin259/rna-velocity-benchmark/resolve/cd1adffe90b778ec1d11ed5985284f1a5d41095c/source_data/v1.0.0/SOURCE_INDEX.tsv?download=true) |
| Array field descriptions | [ARRAY_DICTIONARY.tsv](https://huggingface.co/datasets/Zelin259/rna-velocity-benchmark/resolve/cd1adffe90b778ec1d11ed5985284f1a5d41095c/source_data/v1.0.0/ARRAY_DICTIONARY.tsv?download=true) |

Download only the CSV, TSV or NPZ files you need; individual files require no
ZIP extraction. For all files, use the optional
[Source_Data_complete.zip](https://huggingface.co/datasets/Zelin259/rna-velocity-benchmark/resolve/cd1adffe90b778ec1d11ed5985284f1a5d41095c/source_data/v1.0.0/Source_Data_complete.zip?download=true) (138.15 MB).
It extracts into Source_Data/ with the same directory layout. The
[download manifest](source_data_manifest.tsv) records every individual file
and the optional archive, with immutable URLs, sizes and SHA-256 checksums.

You can also download the directory using the Hugging Face Python client:

```python
from huggingface_hub import snapshot_download

snapshot_download(
    repo_id="Zelin259/rna-velocity-benchmark",
    repo_type="dataset",
    revision="cd1adffe90b778ec1d11ed5985284f1a5d41095c",
    allow_patterns="source_data/v1.0.0/**",
    ignore_patterns="*.zip",
    local_dir="benchmark-results",
)
```

CSV/TSV files retain their saved precision. Load NPZ files with NumPy and
allow_pickle=False; read_arrays.py resolves exact duplicate-array aliases
within each file. The reproduction/fig5 and reproduction/fig6a directories
contain runnable plotting entries and usage instructions. Their inputs and
retained images are available in the same directory tree. The labelled Dynamo
RPE1 entry accepts either the registered H5AD or the included saved
coordinate/velocity arrays.

This correction updates S1 method configurations and DeepKINET dependency
versions, S2 access metadata and the S9 source mapping. Scores, saved plotting arrays and images retain their
published contents. The [original directory snapshot](https://huggingface.co/datasets/Zelin259/rna-velocity-benchmark/tree/f7751b9221683210309d9110d3b55ce842992410/source_data/v1.0.0) remains available. The former size-split ZIP files remain accessible
at their [original fixed commit](https://huggingface.co/datasets/Zelin259/rna-velocity-benchmark/tree/768c153ff25d9e2359d34dc252905e45a7e36cf8/source_data/v1.0.0).
Native method outputs are available through the [result download guide](RESULT_DOWNLOADS.md).
