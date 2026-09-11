# Canonical gene panels

This directory contains the paper-used Gene panels as ordinary repository
files. The panel manifest records each paper-panel mapping and SHA256 value.

All phase-portrait panels use the enlarged canonical display contract:

```text
log1p compression
reference quantile: 0.90
reference length: 0.05
maximum length: 0.10
```

This transformation changes only the display vectors; it does not alter CBDir
or any source H5AD layer. External method results are referenced through logical artifact
locators in `provenance.json`, so the panel files do not depend on a machine
specific filesystem path.
