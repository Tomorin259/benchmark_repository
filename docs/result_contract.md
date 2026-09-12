# Result contract

Method outputs remain method-native H5AD files. A released result must retain
cell and gene identifiers and the layers or embeddings consumed by its
evaluation module. The exact external file, checksum, and required keys are
recorded in `artifacts.tsv`. For DeepKINET, the formal result additionally
records the benchmark seed and deterministic-runtime settings in its method
manifest; those settings are part of the result provenance, not a change to
the upstream model equations.

Evaluation modules write only compact, paper-facing tables and selected panels
to their local `results/` and `figures/` directories. Large H5AD files,
checkpoints, intermediate matrices, and exploratory plots remain external.

Ensemble H5AD artifacts use schema 6. They retain the signed-L1 transition
operator, explicit scVelo UMAP projection, logical input identifiers and
checksums, fixed-neighbor metadata, and aggregation diagnostics. The induced
gene-space `ensemble_velocity = T @ Ms - Ms` exists only to satisfy scVelo's
projection interface; it is not a gene-level aggregation of the four methods.
Build-time PCA from a fixed-neighbor asset is not installed into the result.
