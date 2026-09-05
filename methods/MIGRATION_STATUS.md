# Method-script migration status

This record describes the direct scripts in this directory as of the current
benchmark curation pass. A smoke test is a deliberately small real execution
that writes only into ignored `outputs/` and is removed afterwards; it is not a
substitute for a full scientific rerun.

| Configurations | Direct script | Light validation | Status |
| --- | --- | --- | --- |
| C01–C03 scVelo | `scvelo.py` | `capo_pre` deterministic; compared with the archived workflow | migrated |
| C04 MultiVelo | `multivelo.py` | `capo_pre` RNA-only call | migrated |
| C05–C06 UniTVelo | `unitvelo.py` | `capo_pre` mode 1 | migrated |
| C07 Dynamo | `dynamo.py` | `capo_pre` stochastic Dynamo path | migrated |
| C08 TIVelo | `tivelo.py` | `capo_pre`; verified no redundant second graph/embedding pass | migrated |
| C09 cell2fate | `cell2fate.py` | input contracts checked; model import blocked | script migrated; environment repair pending |
| C10 TFvelo | — | formal output contract and upstream workflow audited | source/package unavailable; do not invent a script |
| C11 VeloAE | `veloae.py` | checkpoint load, forward pass, and projection | migrated |
| C12 VeloVI | `velovi.py` | one-epoch `capo_pre` run, H5AD and stream plot | migrated |
| C13 DeepKINET | `deepkinet.py` | one-epoch `capo_pre` GPU run with raw-count checks; seed-123 deterministic settings recorded | migrated |
| C14 TSvelo | — | formal input/result contract audited | target data contract decision pending |
| C15 cellDancer | `celldancer.py` | three genes / 50 epochs, copyback and graph | migrated |
| C16 DeepVelo (GCN) | `deepvelo.py` | one-epoch `capo_pre` GPU run | migrated; archived env import is broken |
| C17 SDEvelo | `sdevelo.py` | one-epoch `capo_pre` GPU run | migrated; archived env import is broken |
| C18 DeepVelo reconstruction | `deepvelo2022.py` | small real scVelo + TensorFlow one-epoch path | migrated |
| C19 GraphVelo / scVelo | `graphvelo_scv.py` | small real dynamical + tangent-projection path | migrated |
| C20 GraphVelo / Dynamo | `graphvelo_dynamo.py` | small real Dynamo + tangent-projection path | migrated |

## Deliberately unresolved configurations

`cell2fate_env` currently fails before model construction because its pinned
Lightning stack cannot find `wandb`; it also has a Torch version drift from
cell2fate's declared dependency. The direct script is ready, but this shared
environment is not modified during repository curation. Test it after creating
or repairing an isolated compatible environment.

The C10 TFvelo formal artifacts are available, but neither a callable TFvelo
package nor the upstream source checkout is present. Reconstructing an
approximation from saved H5AD fields would not be an official reproduction, so
no substitute script is included.

For C14 TSvelo, the current benchmark `mousebrain.h5ad` contains 3365 cells.
The archived formal `results/TSvelo/mousebrain/TSvelo.h5ad` contains those same
3365 cells plus 288 legacy-only cells (3653 total). Therefore one script cannot
truthfully reproduce both contracts. TSvelo is also only available as a local
source tree, not as an installed package. Before adding C14, choose one of:

1. Preserve the 3653-cell legacy TSvelo artifact and document it as a distinct
   input contract.
2. Rerun TSvelo from the current 3365-cell benchmark object and replace the
   legacy artifact after review.

The second choice is the coherent route for a unified current benchmark, but
it changes a formal result and is intentionally not performed by this script
migration pass.

## Environment notes

The validated execution environments are deliberately not vendored here.
Important deviations discovered during curation are recorded instead:

- `deepvelo_env` lacks `typing_extensions`, preventing PyTorch import. C16's
  smoke test used the compatible `deepvelo2` environment with the same
  DeepVelo/scVelo API family.
- `sdevelo_env` combines Python 3.8 with scVelo 0.3.3, which fails on import.
  C17's smoke test used `sdevelo_py310` with the same SDEvelo 0.2.12 and
  scVelo 0.3.3, with user site packages disabled.
- C18 was tested in `deepvelo2022`, the archived TensorFlow 2.6 / scVelo 0.2.4
  environment. C19 and C20 were tested in `graphvelo`.

No script writes to `/mnt/second19T/zhaozelin/velocity/results`; full reruns
should first write under the ignored `outputs/` area or a separately reviewed
staging location.
