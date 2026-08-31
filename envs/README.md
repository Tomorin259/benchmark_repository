# Environment evidence

Supplementary Table S1 reports the package and selected dependency versions used in the benchmark. `manifests/environments.tsv` connects those statements to Conda environments still present on the audited server.

Most assignments are inferred from matching environment names and package tuples. They are not lock files. Only DeepKINET, TSvelo, and DeepVelo_scVelo have nearby method-specific specifications, and each has a material limitation: version drift, a missing source revision, conflicting installed metadata, or unresolved asset licensing.

For that reason, Phase 2 does not export the current Conda environments. A raw `conda env export` would capture unrelated build strings, editable installs, local paths, and historical contamination without proving which packages produced the paper results.

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
