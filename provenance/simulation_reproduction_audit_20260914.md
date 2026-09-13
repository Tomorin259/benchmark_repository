# Simulation reproduction audit

This note records the audit of the three generated simulation inputs used by
the gene-level benchmark. It is an audit record, not a generator: the private
source scripts and temporary reproduction outputs are not added to this
repository.

## Audited source and scenarios

The audited source consisted of the following scripts (the paths are the
repository-relative names used for source identification; these scripts are
not currently part of the public checkout):

| Scenario | Audited script | Canonical input |
| --- | --- | --- |
| murk / one-sided boost | `indicators/gene/simulation/test_only_murk.py` | `dataset-sim-only-murk` |
| branching | `indicators/gene/simulation/test_only_branching.py` | `dataset-sim-only-branching` |
| hybrid/all | `indicators/gene/simulation/test_all.py` | `dataset-sim-hybrid` |

The workflow uses `cellDancer` 1.1.7 through its public Python API
(`celldancer.simulation`); cellDancer source is not vendored. Its dependency
license and the exact upstream source revision still need to be recorded in a
future release environment/third-party notice.

## RNG, equations, and deterministic boundary

The audited scripts do not set an explicit seed. In the audited cellDancer
implementation, each simulated gene consumes NumPy normal draws in
`_jitter`, while the scripts also consume uniform draws for the one-sided
boost parameters. Adding an explicit `np.random.seed(seed)` before generation
makes a given installed implementation deterministic: repeated generation
with the same seed reproduced the same arrays, IDs, and metadata.

The generated truth equations were checked directly:

```text
velocity_s_true = beta_true * unspliced - gamma_true * spliced
velocity_u_true = alpha_true - beta_true * unspliced
```

Both residuals were zero in the audited outputs. This establishes structural
and numerical consistency of the generation equations, not identity with the
historical canonical files.

## Canonical inputs

The following are the canonical files already represented by the corresponding
rows in the root `artifacts.tsv` registry:

| Logical artifact | Shape | Size (bytes) | SHA256 |
| --- | ---: | ---: | --- |
| `dataset-sim-only-murk` | 1000 x 650 | 39,807,648 | `e853780b3cb28ae9941972950e0b02098f57b5f7e35908ab0d5d2473fe0e6bc3` |
| `dataset-sim-only-branching` | 2000 x 620 | 75,781,536 | `392eccf96cd31f27e8b99327f0b061e475d78a1cfea54c75af4564ddb163fe14` |
| `dataset-sim-hybrid` | 2000 x 720 | 87,805,952 | `9b9bebf9159e02c1bcb880d8aba323982e70b7faa67d6d92c6313db603e72ee2` |

The public registry already contains these hashes; no registry update was
needed for this audit.

## Reproduction checks

Full isolated workflows with candidate seeds 0, 42, and 123 were run for all
three scenarios. None matched the canonical numerical arrays (and therefore
none was byte-exact). The generated files retained the canonical structural
shape, observation/variable IDs, observation metadata, `X_pca`/`X_umap`, and
the expected `uns` keys, but their numerical arrays differed.

For the murk scenario, the first canonical one-sided-boost parameters were
also compared after consuming the expected housekeeping normal draws. A
search of startup seeds 0 through 1000 produced no match, so the historical
RNG state cannot be explained by a simple seed in that tested range.

For reference, the isolated seed-123 output hashes were:

| Scenario | Audit-only output SHA256 |
| --- | --- |
| murk | `394e336c438358c7a7d176f1142ab2eb5820eb90e874e11c26a869f6014db91c` |
| branching | `f8721a0030b71e179ed21cfcf13ed0a6873b5e40ac76433abeb36b3b1699ffd5` |
| hybrid/all | `e51000428f75c3f18018fde155bafb8e5dab15e13533c7b6797ef02e94ede294` |

These hashes are temporary audit evidence only and are not canonical artifact
hashes or published inputs.

## Known limitations and source bugs

- No invocation record, explicit seed, environment lock, or upstream commit is
  present in the original scripts.
- `test_only_branching.py` and `test_all.py` both write an output name inherited
  from the murk workflow, `s1_only_murk.h5ad`, despite representing the
  branching and hybrid/all scenarios.
- The scripts combine simulation, dimensionality reduction, preprocessing,
  and plotting, so they do not provide a clean generation entry point.
- `test_all.py` sets `n_multi_bwd = 0`; the nominal MultiBackward component is
  therefore not materialized by that workflow.
- The audited scripts depend on external source and local path conventions;
  those paths are not a public reproducibility contract.

## Conclusion

The audit supports a claim of structural reproduction and equation-level
consistency under cellDancer 1.1.7 with an explicitly controlled RNG. It does
not support a claim that the historical canonical H5AD files can be regenerated
byte-for-byte, or that their exact initial RNG state is recoverable from the
available scripts. Public documentation should therefore describe these
inputs as audited generated artifacts, not as byte-exact reproducible outputs.
