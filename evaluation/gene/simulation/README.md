# Simulation input generation

`generate.py` is the public generation entry point for the three simulation
scenarios used by the gene-level benchmark. It contains simulation only;
preprocessing, dimensionality reduction, and plotting remain separate
downstream operations.

To fit methods on the three saved benchmark inputs, follow the nine
[simulation method commands](../../../methods/README.md#run-methods-on-simulations).
The [Gene evaluation guide](../README.md#recalculate-cbdir) then recalculates
CBDir and CosSim from the method outputs.

The scientific parameters live in [`config.json`](config.json), not in the
command line. The entries preserve the audited CellDancer 1.1.7 logic for
`murk`, `branching`, and `hybrid`. The branching BoostALL parameters are
`alpha1=5, alpha2=4, beta1=1, beta2=2, gamma1=5, gamma2=1`.
These values follow the canonical H5AD parameter layers. Simulator source
hashes and the parameter-source record are retained in `config.json`.

The generator requires cellDancer 1.1.7 and `setuptools<81`, which provides
`pkg_resources`. Run from the repository root with a relative output directory:

```bash
python evaluation/gene/simulation/generate.py \
  --scenario branching --output-root outputs/simulation
```

The output has the correct scenario-specific name (`sim_only_murk.h5ad`,
`sim_only_branching.h5ad`, or `sim_hybrid.h5ad`) and a sidecar
`*.invocation.json` containing its config hash, generator hash, output hash,
shape, dependency versions, random-stream policy, and canonical artifact ID.
The benchmark uses the historical simulation H5ADs indexed in
[`artifacts.tsv`](../../../artifacts.tsv). This generator creates additional
simulated datasets. Its randomness settings are recorded in `config.json`
and each output's invocation record.
