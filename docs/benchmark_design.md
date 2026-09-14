# Benchmark design

The benchmark is organized by the scientific questions reported in the paper:

1. cell-level trajectory inference across 14 real datasets;
2. gene-level velocity evaluation on real and simulated data;
3. kinetic-parameter recovery on RPE1 and organoid scEU-seq data; and
4. post-hoc signed-transition ensemble evaluation.

Method reproduction and evaluation are separate stages. Scripts in
[`../methods/`](../methods/) produce method-native results. Each module under
[`../evaluation/`](../evaluation/) owns its metric definition, input mapping,
alignment rules, compact results, and scientific panels.

The ensemble consumes released method results after method execution. It is
not an additional registered method. Its comparator definitions and evaluation
contract are owned by [`../ensemble/`](../ensemble/).
