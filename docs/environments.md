# Environment evidence and release requirements

Supplementary Table S1 reports the package and selected dependency versions used in the benchmark. The repository currently contains no reviewed Conda specification, lock file, or container definition. This document records the requirements that future reviewed environment files must satisfy. Historical environment assignments and import checks remain in the private audit records.

Most assignments are inferred from matching environment names and package tuples. They are not lock files. The Phase 2 audit found an importable candidate for each of the 16 configuration environments; these checks were package or entry-point import smoke tests and did not rerun model training. DeepVelo GCN imports with a NumPy/SciPy compatibility warning. The SDEvelo environment used for the original check is not import-compatible, so the audited C17 candidate is the working `sdevelo_py310` environment; its release specification still requires curation. DeepKINET, veloAE, and DeepVelo_scv also require specification cleanup because of version drift, editable or historical metadata, or missing revision information.

For that reason, Phase 2 does not export the current Conda environments. A raw `conda env export` would capture unrelated build strings, editable installs, local paths, and historical contamination without proving which packages produced the paper results. A passing import smoke test must not be read as a complete reproduction claim.

Historical evidence, runtime smoke results, and public-release readiness are
intentionally separate facts. They belong to the private audit record until a
clean specification is reviewed and published.

## Requirements for a reviewed specification

Each future environment file must:

- identify exactly one or an explicitly shared set of configuration IDs;
- pin the benchmarked method version or upstream commit;
- reproduce the manuscript-reported Python and key dependency versions, or document a validated deviation;
- avoid local filesystem paths and unreviewed editable installs;
- record GPU/CUDA constraints where required;
- install successfully in a clean environment;
- pass the corresponding runner smoke test;
- preserve third-party license notices without vendoring unapproved assets.

The audited Conda paths are provenance evidence only and will not appear in the final public execution interface.
