# Naming contract

This file is the single source of truth for names used in manuscript-facing
figures, captions, tables, and compact plot labels. Operational IDs and result
paths remain stable; display labels are applied only at the presentation
boundary.

## Method labels

Use this order and spelling for the 20 benchmark configurations:

| # | Operational ID | Manuscript/display label |
| ---: | --- | --- |
| 1 | `scvelo_deterministic` | `scVelo_deterministic` |
| 2 | `scvelo_stochastic` | `scVelo_stochastic` |
| 3 | `scvelo_dynamic` | `scVelo_dynamic` |
| 4 | `multivelo` | `MultiVelo (RNA_only)` |
| 5 | `unitvelo` | `UniTVelo` |
| 6 | `unitvelo_mode2` | `UniTVelo(mode2)` |
| 7 | `dynamo` | `Dynamo(RNA_seq)` |
| 8 | `TIvelo` | `TIvelo` |
| 9 | `cell2fate` | `cell2fate` |
| 10 | `TFvelo` | `TFvelo` |
| 11 | `veloae` | `veloAE` |
| 12 | `veloVI` | `veloVI` |
| 13 | `DeepKINET` | `DeepKINET` |
| 14 | `TSvelo` | `TSvelo` |
| 15 | `celldancer` | `cellDancer` |
| 16 | `deepvelo` | `DeepVelo` |
| 17 | `SDEvelo` | `SDEvelo` |
| 18 | `deepvelo2022` | `DeepVelo_scv` |
| 19 | `graphvelo_scv` | `GraphVelo_scv` |
| 20 | `graphvelo_dynamo` | `GraphVelo_dyn` |

Display rules:

- Do not append `Velocyto` or `Velocyto model` to `scVelo_deterministic`.
- Keep the `(RNA_only)` qualifier in `MultiVelo (RNA_only)`.
- Keep `DeepVelo`, `DeepVelo_scv`, `GraphVelo_scv`, and `GraphVelo_dyn`
  distinct.
- Preserve capitalization and punctuation exactly as listed.
- Do not rename operational IDs, result directories, or source filenames.

## Real-dataset labels

Use these display labels for the 14 real datasets:

| Internal ID or evaluation alias | Manuscript/display label |
| --- | --- |
| `pancreas` | `Pancreas` |
| `gastrulation` | `Gastrulation` |
| `forebrain` | `ForebrainGlut` |
| `DentateGyrus_lamanno`, `dentategyrus_lamanno` | `Dentategyrus_lamanno` |
| `dentategyrus` | `Dentategyrus` |
| `hematopoiesis` | `Human hematopoiesis` |
| `bonemarrow` | `Human bone marrow` |
| `rpe1`, `sceu_rpe1` | `rpe1` |
| `organoid`, `sceu_organoid` | `organoid` |
| `capofibro`, `capo_pre` | `capofibro (VeloCycle)` |
| `rpe1_gen` | `rpe1 (VeloCycle)` |
| `retina`, `Mouse_retina` | `Mouse retina` |
| `mousebrain` | `mousebrain` |
| `neuron`, `scnt_neuron` | `scNT_neuron` |

`rpe1` and `rpe1 (VeloCycle)` are different datasets and must never be merged
or treated as display aliases. The first is the Battich scEU-seq RPE1 dataset;
the second is the Lederer/VeloCycle RPE1 dataset. `capo_pre` and `capofibro`
are internal aliases for `capofibro (VeloCycle)`.

## Simulated-dataset labels

| Scenario | Manuscript/display label |
| --- | --- |
| normal genes plus MURK genes | `Normal + MURK` |
| normal genes plus branching genes | `Normal + branching` |
| normal, MURK, and branching genes | `Hybrid (normal + MURK + branching)` |

## Use in code and artifacts

- Metric CSVs, manifests, configuration files, result paths, and source
  filenames use operational IDs.
- Figures, captions, manuscript tables, and prose use display labels.
- Any display-name mapping in plotting code must follow this contract exactly; do not introduce a conflicting
  spelling table.
- Historical snapshots may retain their original labels; all current and new
  artifacts must follow this contract.
