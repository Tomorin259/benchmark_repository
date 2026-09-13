# Figure map

Figures are generated beside the experiment that owns their data and metric
definition. Each experiment keeps its panel assets and reproducible plot code
inside its own module; the repository-level figures/ directory contains only
confirmed paper-level composite PDFs:

- cell-level panels come from `evaluation/cell/`;
- gene-level real and simulation panels come from `evaluation/gene/`;
- kinetics-recovery panels come from `evaluation/kinetics_recovery/plotting/`; and
- integrated summaries such as Figure 6 come from `evaluation/summary/`.

For Figure 2, cell-level panels are in
`evaluation/cell/current/figures/`. For Figure 5, the panel assets are in
`evaluation/kinetics_recovery/figures/fig5/`, with
`plotting/figure5.py` owning the legacy-style quantitative panels from the
explicit `plotting/figure5_config.json` input map and
`plotting/streams.py` owning organoid streams and the current DeepKINET RPE1
stream. The other four RPE1 stream images in Figure 5A remain recorded manual
assets. The input map keeps the historical cellDancer/DeepVelo/Dynamo/TSvelo
data results and switches only the DeepKINET gamma results to the final
seed-123 generation. The two kinetic-class Venn panels are recomputed from
their mapped alpha/gamma trajectories and write per-gene membership tables.
The original paper-level
PDFs are in `figures/`. When a manuscript figure combines panels from
multiple modules, this file records the assembly order and source panel names.

The assembly manifest records five known panel-level source gaps for follow-up:
the Figure 1 study-design panels A and B are missing, the Figure 3 simulation
schematic is missing, the Figure 4 branching overview source is unverified, and
the Figure 6 integrated panel source is missing. Supplementary Figure 9 now
contains reproducible RowGeometricMedian-4 CBDir and BATC metric panels. These
are documented post-release source-recovery items; no private path or
replacement source is inferred.
