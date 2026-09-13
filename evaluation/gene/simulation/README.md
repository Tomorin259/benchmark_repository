# Simulation input generation

`generate.py` is the public generation entry point for the three simulation
scenarios used by the gene-level benchmark. It contains simulation only;
preprocessing, dimensionality reduction, and plotting remain separate
downstream operations.

The scientific parameters live in [`config.json`](config.json), not in the
command line. The entries preserve the audited CellDancer 1.1.7 logic for
`murk`, `branching`, and `hybrid`. The branching BoostALL parameters are
`alpha1=5, alpha2=4, beta1=1, beta2=2, gamma1=5, gamma2=1`.
The archived historical branching script contained an older Branch-B setting;
the values above come from the canonical H5AD parameter layers and were
validated by a full regenerated branching workflow. Hashes of all three
audited source scripts are retained in `config.json`.

Run from the repository root with a relative output directory:

```bash
python evaluation/gene/simulation/generate.py \
  --scenario branching --output-root outputs/simulation
```

The output has the correct scenario-specific name (`sim_only_murk.h5ad`,
`sim_only_branching.h5ad`, or `sim_hybrid.h5ad`) and a sidecar
`*.invocation.json` containing its config hash, generator hash, output hash,
shape, dependency versions, random-stream policy, and canonical artifact ID.
No seed is set: random kinetic parameters and CellDancer simulation noise are
part of the study design, so every invocation is an independent realization.
Consequently a new output is not expected to be byte-identical to either a
second run or the historical canonical H5AD.

The historical inputs remain the artifacts indexed by the root
[`artifacts.tsv`](../../../artifacts.tsv). Generated outputs from this entry
point are new audit/regeneration products and must not be registered as
canonical replacements until independently validated and published.
