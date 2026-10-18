# Names and identifiers

| Registry | Identifier fields | Display field |
| --- | --- | --- |
| [Method configurations](../methods/configurations.tsv) | `config_id` (C01–C20), `operational_id` | `manuscript_name` |
| [Datasets](../data/datasets.tsv) | Stable dataset ID | `paper_name` |

Use stable IDs for paths and joins, and registry display fields in figures,
tables and the website. Evaluation-specific aliases are resolved at input
boundaries.

DeepVelo and DeepVelo_scv are distinct configurations; GraphVelo_scv and
GraphVelo_dyn use different upstream velocities. The operational ID `velovi`
has the display spelling **veloVI**.
The [Mean-4 and GeoMedian-4 ensembles](../ensemble/README.md) are post-hoc
comparators outside the 20-configuration registry.
