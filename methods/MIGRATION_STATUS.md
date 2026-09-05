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
| C10 TFvelo | `tfvelo.py` | official source import/help; 14 existing result contracts | source-backed launcher; full rerun pending |
| C11 VeloAE | `veloae.py` | checkpoint load, forward pass, and projection | migrated |
| C12 VeloVI | `velovi.py` | one-epoch `capo_pre` run, H5AD and stream plot | migrated |
| C13 DeepKINET | `deepkinet.py` | one-epoch `capo_pre` GPU run with raw-count checks; seed-123 deterministic settings recorded | migrated |
| C14 TSvelo | `tsvelo.py` | official source import/help; 14 existing result contracts | source-backed launcher; full rerun pending |
| C15 cellDancer | `celldancer.py` | three genes / 50 epochs, copyback and graph | migrated |
| C16 DeepVelo (GCN) | `deepvelo.py` | one-epoch `capo_pre` GPU run | migrated; archived env import is broken |
| C17 SDEvelo | `sdevelo.py` | one-epoch `capo_pre` GPU run | migrated; archived env import is broken |
| C18 DeepVelo reconstruction | `deepvelo2022.py` | small real scVelo + TensorFlow one-epoch path | migrated |
| C19 GraphVelo / scVelo | `graphvelo_scv.py` | small real dynamical + tangent-projection path | migrated |
| C20 GraphVelo / Dynamo | `graphvelo_dynamo.py` | small real Dynamo + tangent-projection path | migrated |

## Source-backed caveats

`cell2fate_env` currently fails before model construction because its pinned
Lightning stack cannot find `wandb`; it also has a Torch version drift from
cell2fate's declared dependency. The direct script is ready, but this shared
environment is not modified during repository curation. Test it after creating
or repairing an isolated compatible environment.

The C10 TFvelo source checkout is available at
`/mnt/second19T/zhaozelin/TFvelo-main`. Its official demo scripts use relative
database paths and the analysis entry point hard-codes the hematopoiesis
`rc.h5ad`; `methods/tfvelo.py` stages each configured input and patches only
that path in a temporary analysis copy. The model implementation remains in
the official checkout.

The C14 TSvelo source checkout is available at
`/mnt/second19T/zhaozelin/velocity/methods/TSvelo`. Its official preprocessing
falls back to an organoid path for unrecognized dataset names;
`methods/tsvelo.py` stages each configured input at that compatibility path
and invokes the original `TSvelo_run.py` entry point.

The TSvelo mousebrain result intentionally uses the same biological dataset
through a different observation panel: the 4881-cell multiome loom is reduced
to 3653 cells by the official preprocessing. Cell-level evaluation maps those
cells to the common 3365-cell benchmark panel by cell ID; the method-native
3653-cell H5AD is preserved unchanged.

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
- C10 uses the external TFvelo source checkout with the audited
  `TFvelo_env` versions (scVelo 0.2.4, Scanpy 1.9.3, AnnData 0.8.0,
  NumPy 1.21.6, SciPy 1.10.1); TFvelo is imported from source rather than an
  installed distribution.
- C14 uses `TSvelo_env` with the pinned versions in its environment file,
  including scVelo 0.3.2, Scanpy 1.9.8, AnnData 0.9.2, and PyTorch 2.4.1.

No script writes to `/mnt/second19T/zhaozelin/velocity/results`; full reruns
should first write under the ignored `outputs/` area or a separately reviewed
staging location.
