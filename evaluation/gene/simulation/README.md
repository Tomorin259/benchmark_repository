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
  --scenario branching --seed 123 --output-root outputs/simulation
```

The output has the correct scenario-specific name (`sim_only_murk.h5ad`,
`sim_only_branching.h5ad`, or `sim_hybrid.h5ad`) and a sidecar
`*.invocation.json` containing the seed, config hash, shape, command, and
canonical artifact ID. The seed makes a new run deterministic under the
reviewed NumPy/CellDancer versions; it does not make the historical canonical
H5AD byte-exactly reproducible.

The audited environment generated each scenario twice with seed 123 and
obtained byte-identical H5AD files. The expected hashes for that regression
run are stored in `config.json`; every invocation sidecar records the actual
output, configuration, generator, and dependency-version hashes or values.

The historical inputs remain the artifacts indexed by the root
[`artifacts.tsv`](../../../artifacts.tsv). Generated outputs from this entry
point are new audit/regeneration products and must not be registered as
canonical replacements until independently validated and published.
