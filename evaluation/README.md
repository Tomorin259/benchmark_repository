# Evaluation

Each module provides its input map, metric implementation, compact results and
plotting entry points.

| Task | Results and commands |
| --- | --- |
| Cell-level velocity evaluation on 14 real datasets | [Cell](cell/README.md) |
| Gene-level velocity directions on real and simulated data | [Gene](gene/README.md) |
| Rates against RPE1 and organoid scEU-seq references | [Kinetic-parameter recovery](kinetics_recovery/README.md) |
| Figure 7 and cross-experiment summaries | [Summary](summary/README.md) |
| Cell-level comparison on the three simulations | [Simulation Cell](cell/simulation/README.md) |
| Post-hoc signed-transition aggregation | [Ensemble](../ensemble/README.md) |

Start with committed tables for inspection, or materialize the method-native
H5ADs and reference inputs to recalculate metrics. See
[data access](../data/README.md) and [environments](../docs/environments.md).
