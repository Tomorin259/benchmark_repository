# Result contract

Method execution produces method-native H5AD files. A released result must
retain cell and gene identifiers and every layer, embedding, graph, or metadata
field consumed by its evaluation module.

Result metadata has separate authoritative locations:

- [`../artifacts.tsv`](../artifacts.tsv) records the logical artifact ID,
  tracked path when applicable, SHA256, release status, stable URI, and notes;
- [`../methods/configurations.tsv`](../methods/configurations.tsv) records the
  registered method configuration and its standard result layer and embedding;
- evaluation input tables record the exact result file, keys, analysis role,
  and alignment information used by that evaluation;
- method run metadata records execution-specific settings such as a fixed seed
  or deterministic runtime mode.

Evaluation modules write compact, paper-facing tables, audit records, and
selected panels to their local `results/` and `figures/` directories. Large
H5AD files, checkpoints, intermediate matrices, and exploratory plots remain
external artifacts.

Specialized result schemas belong to their producing module. The signed
transition, fixed-neighbor, projection, and diagnostic fields of ensemble H5AD
files are defined in [`../ensemble/README.md`](../ensemble/README.md).
