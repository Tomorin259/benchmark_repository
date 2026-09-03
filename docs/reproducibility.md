# Reproducibility policy

Method scripts follow the official tutorial and reviewed environment for that
method. The repository does not add a global seed unless the upstream workflow
specifies one. Environment names, local editable checkouts, server paths, and
historical result generations are not treated as reproducibility contracts.

When a method requires raw counts, its own script documents the selected
`_counts.h5ad` artifact, alignment, count checks, and any preprocessing rule.
This requirement is local to the consuming method.
