# Naming conventions

Machine-readable registries are the source of truth for names. Do not copy
their complete contents into prose or infer public names from Python filenames.

## Methods

[`../methods/configurations.tsv`](../methods/configurations.tsv) is the sole
registry of the 20 benchmark configurations:

- `config_id` is the stable C01--C20 record identifier;
- `operational_id` identifies the runnable configuration and result family;
- `manuscript_name` is the display name used in paper tables, figures, and
  prose.

Operational IDs are used for paths and joins. Presentation code must obtain
display names from `manuscript_name` rather than maintaining a private rename
dictionary.

Important distinctions retained by the registry include:

- `velovi` is the operational ID, while **veloVI** is the display spelling;
- C16 DeepVelo and C18 DeepVelo_scv are distinct configurations;
- C19 GraphVelo_scv and C20 GraphVelo_dyn differ by their upstream velocity
  input.

## Datasets

[`../data/datasets.tsv`](../data/datasets.tsv) owns stable dataset IDs and
paper-facing display names. Evaluation-specific aliases are resolved at the
evaluation input boundary; they must not replace the registered dataset ID.

## Ensemble comparators

Signed-T Mean-4 and Signed-T RowGeometricMedian-4 are post-hoc comparators, not
method registry rows. Their identifiers, definitions, and ranking rules are
owned by [`../ensemble/`](../ensemble/).
