# Figure map

Figures are generated beside the experiment that owns their data and metric
definition:

- cell-level panels come from `evaluation/cell/`;
- gene-level real and simulation panels come from `evaluation/gene/`;
- kinetics-recovery panels come from `evaluation/kinetics_recovery/`; and
- integrated summaries such as Figure 6 come from `evaluation/summary/`.

When a manuscript figure combines panels from multiple modules, this file
records the assembly order and source panel names. The final manually composed
paper PDF is not required to be stored in the code repository.
