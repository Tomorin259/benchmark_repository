# Result contract

Method outputs remain method-native H5AD files. A released result must retain
cell and gene identifiers and the layers or embeddings consumed by its
evaluation module. The exact external file, checksum, and required keys are
recorded in `artifacts.tsv`.

Evaluation modules write only compact, paper-facing tables and selected panels
to their local `results/` and `figures/` directories. Large H5AD files,
checkpoints, intermediate matrices, and exploratory plots remain external.
