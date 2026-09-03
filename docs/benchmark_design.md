# Benchmark design

The benchmark is organized by the scientific question reported in the paper:

1. cell-level trajectory inference across 14 real datasets;
2. gene-level velocity evaluation on real and simulated data;
3. kinetic-parameter recovery on RPE1 and organoid scEU-seq data; and
4. an optional post-hoc signed-transition ensemble.

Method reproduction is kept separate from evaluation. Each evaluation module
reads the native result files produced by the methods and translates their
method-specific keys at the input boundary.
