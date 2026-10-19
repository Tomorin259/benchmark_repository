# Source Data

[Browse Source Data v1.0.0](https://huggingface.co/datasets/Zelin259/rna-velocity-benchmark-results/tree/768c153ff25d9e2359d34dc252905e45a7e36cf8/source_data/v1.0.0)
for the main and supplementary figures and tables. The seven ZIP files contain
88 numerical/metadata tables, 243 saved plot-input arrays, plotting scripts,
and retained image materials.

| Download | Contents | Size |
| --- | --- | ---: |
| [Source_Data.zip](https://huggingface.co/datasets/Zelin259/rna-velocity-benchmark-results/resolve/768c153ff25d9e2359d34dc252905e45a7e36cf8/source_data/v1.0.0/Source_Data.zip?download=true) | Tables, indexes and plotting scripts | 1.47 MB |
| [Source_Data_display_01.zip](https://huggingface.co/datasets/Zelin259/rna-velocity-benchmark-results/resolve/768c153ff25d9e2359d34dc252905e45a7e36cf8/source_data/v1.0.0/Source_Data_display_01.zip?download=true) | Saved display arrays | 25.25 MB |
| [Source_Data_display_02.zip](https://huggingface.co/datasets/Zelin259/rna-velocity-benchmark-results/resolve/768c153ff25d9e2359d34dc252905e45a7e36cf8/source_data/v1.0.0/Source_Data_display_02.zip?download=true) | Saved display arrays | 25.61 MB |
| [Source_Data_display_03.zip](https://huggingface.co/datasets/Zelin259/rna-velocity-benchmark-results/resolve/768c153ff25d9e2359d34dc252905e45a7e36cf8/source_data/v1.0.0/Source_Data_display_03.zip?download=true) | Saved display arrays | 26.13 MB |
| [Source_Data_display_04.zip](https://huggingface.co/datasets/Zelin259/rna-velocity-benchmark-results/resolve/768c153ff25d9e2359d34dc252905e45a7e36cf8/source_data/v1.0.0/Source_Data_display_04.zip?download=true) | Saved display arrays | 26.61 MB |
| [Source_Data_display_05.zip](https://huggingface.co/datasets/Zelin259/rna-velocity-benchmark-results/resolve/768c153ff25d9e2359d34dc252905e45a7e36cf8/source_data/v1.0.0/Source_Data_display_05.zip?download=true) | Saved display arrays | 8.85 MB |
| [Source_Data_materials.zip](https://huggingface.co/datasets/Zelin259/rna-velocity-benchmark-results/resolve/768c153ff25d9e2359d34dc252905e45a7e36cf8/source_data/v1.0.0/Source_Data_materials.zip?download=true) | Original images and reproduction previews | 24.20 MB |

Unzip Source_Data.zip and the display parts you need into the same parent
folder. All parts merge into Source_Data/. Start with SOURCE_INDEX.tsv;
PACKAGE_INDEX.tsv identifies the archive containing each file. The materials
archive contains original images and checked reproduction previews.

CSV/TSV files retain their saved precision. Load NPZ files with NumPy and
allow_pickle=False; read_arrays.py resolves exact duplicate-array aliases.
The reproduction/fig5 and reproduction/fig6a directories contain runnable
entries and usage instructions. The labelled Dynamo RPE1 entry accepts either
the registered H5AD or the included saved coordinate/velocity arrays.

[Download sizes and SHA-256 checksums](source_data_manifest.tsv) identify every
archive. Links are pinned to the published content commit. Native method
outputs remain available through the [result download guide](RESULT_DOWNLOADS.md).
