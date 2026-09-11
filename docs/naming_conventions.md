# Naming conventions

methods/configurations.tsv is the sole 20-method benchmark registry. It
contains exactly 20 unique configuration records. The columns have separate
responsibilities:

- config_id is the stable C01-C20 record identifier.
- operational_id identifies the runnable configuration and result family.
- manuscript_name is the public display name used in tables, figures, and
  captions. The display spelling is veloVI, not velovi.

## Registered configurations

| config_id | operational_id | display name |
| --- | --- | --- |
| C01 | scvelo_deterministic | scVelo_deterministic |
| C02 | scvelo_stochastic | scVelo_stochastic |
| C03 | scvelo_dynamic | scVelo_dynamic |
| C04 | multivelo | MultiVelo (RNA_only) |
| C05 | unitvelo | UniTVelo |
| C06 | unitvelo_mode2 | UniTVelo(mode2) |
| C07 | dynamo | Dynamo(RNA_seq) |
| C08 | TIvelo | TIvelo |
| C09 | cell2fate | cell2fate |
| C10 | TFvelo | TFvelo |
| C11 | veloae | veloAE |
| C12 | velovi | veloVI |
| C13 | deepkinet | DeepKINET |
| C14 | tsvelo | TSvelo |
| C15 | celldancer | cellDancer |
| C16 | deepvelo | DeepVelo |
| C17 | SDEvelo | SDEvelo |
| C18 | deepvelo2022 | DeepVelo_scv |
| C19 | graphvelo_scv | GraphVelo_scv |
| C20 | graphvelo_dynamo | GraphVelo_dyn |

Names and settings follow the existing registry and paper/source evidence;
they are not inferred from Python filenames. In particular, velovi.py is the
implementation filename and is not changed to alter the public veloVI display
spelling.

DeepVelo C16 is the graph-convolutional configuration. C18 is the distinct
DeepVelo reconstruction configuration. GraphVelo C19 consumes pre-calculated
scVelo input, while C20 consumes pre-calculated Dynamo input; these are two
distinct configurations of the same software.

Dataset IDs identify biological inputs and their evaluation panels. They are
not method IDs, configuration IDs, or ensemble IDs.

Mean-4 and CenterWeighted-4 are independent comparators, not registry rows.
Mean-4 is the baseline; CenterWeighted-4 is the formal ensemble. Each is
ranked separately against the 20 registered methods and they are not compared
with each other. Both implementations are in ensemble/. The shared
evaluation/cell/reference_neighbors/ fixed-neighbor contract is used by the
ensemble evaluation; formal fixed-neighbor all-14 results are currently
pending.
