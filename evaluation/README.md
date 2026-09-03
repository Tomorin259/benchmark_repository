# Benchmark evaluation

Each submodule owns one scientific question and keeps its extraction,
dataset-specific processing, metrics, compact results, and generated panels
together. Evaluation consumes method-native H5AD files through explicit key
mappings; it does not impose a universal result schema on the methods.

- `cell/`: cell-level trajectory metrics across the 14 real datasets.
- `gene/`: gene-level metrics for real and simulated data.
- `parameter_recovery/`: RPE1 and organoid kinetic-parameter recovery.
- `summary/`: cross-experiment tables and integrated panels such as Figure 6.
