# Reproducibility policy

Method scripts follow the official tutorial and reviewed environment for that
method. The repository does not add a global seed by default. The formal
DeepKINET benchmark result is the explicit exception: seed 123 and deterministic
PyTorch/CUDA settings are recorded because the upstream workflow leaves the
global RNG state unspecified and the result otherwise depends on run order.
Environment names, local editable checkouts, server paths, and historical
result generations are not treated as reproducibility contracts.

When a method requires raw counts, its own script documents the selected
`_counts.h5ad` artifact, alignment, count checks, and any preprocessing rule.
This requirement is local to the consuming method.

The ensemble pins scVelo 0.3.4 and uses its tutorial graph call with
`vkey="velocity"` and `xkey="Ms"`. Its signed transition call explicitly uses
`scale=10`, self transitions, negative cosines, zero diffusion weight, and the
default transition-neighbor count. The repository-local fixed-neighbor NPZ
files are the public graph artifacts: their numerical arrays match the audited
upstream assets, while their metadata is deliberately path-free and therefore
their file SHA256 values are benchmark-specific.

The ensemble result metadata keeps provenance dimensions separate:
`lineage/current_upstream` identifies `transition_ensemble/main@4f2a621`,
`release_index_commit` identifies the upstream release index commit, and
`numerical_reproduction_commit` identifies the clean source used for numerical
reproduction. These fields must not be collapsed into one `source_commit`.
The detailed evidence is recorded in
[`provenance/fixed_neighbor_reproduction_20260912.md`](../provenance/fixed_neighbor_reproduction_20260912.md).
