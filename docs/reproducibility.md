# Reproducibility policy

Reproduction claims must identify the scientific inputs, method configuration,
software source, randomness policy, output, and validation performed. A
structural rerun must not be described as byte-exact reproduction.

## Method runs

Method scripts follow their audited upstream workflows. Versions, source
references, scientific options, result layers, and embeddings are recorded in
[`../methods/configurations.tsv`](../methods/configurations.tsv). External
source checkouts must use a recorded release or commit before public release.

The repository does not impose one benchmark-wide random seed. Each launcher
must state whether randomness is fixed or intentionally inherited from the
upstream workflow. The formal DeepKINET configuration is the explicit
seed-123, deterministic PyTorch/CUDA case. Raw-count selection and preprocessing
rules remain local to the methods that consume those inputs.

## Data and artifacts

Large inputs and method results may remain outside Git. Their logical IDs,
checksums, publication status, and stable URIs are recorded in
[`../artifacts.tsv`](../artifacts.tsv). Public code must resolve an artifact or
user-supplied root and must not depend on a private server path.

The archived simulation workflow is intentionally stochastic. Its scripts,
canonical output hashes, structural checks, and non-byte-exact boundary are
recorded in
[`simulation_reproduction_audit_20260914.md`](simulation_reproduction_audit_20260914.md).

## Specialized workflows

Module-specific numerical settings belong with the module that executes them.
The ensemble implementation contract is documented in
[`../ensemble/README.md`](../ensemble/README.md). Durable evidence separating
its upstream lineage, release index, and numerical reproduction is recorded in
[`fixed_neighbor_reproduction_20260912.md`](fixed_neighbor_reproduction_20260912.md).

Kinetics reproduction uses the frozen central-result inputs and runtime
versions recorded in `evaluation/kinetics_recovery/results/central_results/manifest.json`.
DeepKINET seed identity is backed by H5AD metadata and input checksums; CSV
byte identity is not implied; small full-table differences have an unestablished cause. Figure 6 streams reuse
saved projections, and quantitative plots consume one evaluator output.
