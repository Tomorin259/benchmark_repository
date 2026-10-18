# Obtain the original kinetic-reference tables

The input distribution provides 25 mirrored scientific files. Obtain the three
additional kinetic-reference dependencies, Battich Tables S1/S3/S4, from the
original study using the steps below.

1. Visit https://doi.org/10.1126/science.aax3072 and open Supplementary Materials.
2. Obtain Tables S1, S3 and S4 through the original publisher/author access route.
3. Place the benchmark-compatible CSV files at the paths in
   `upstream_references.tsv`, using the layout required by your evaluator.
4. From this repository root, run `python data/check_references.py --artifact-root /path/to/artifacts`.
   In the HF version folder, run `python check_references.py --artifact-root /path/to/artifacts` instead.

There are two documented layouts: velocity-root paths used by existing
launchers and `references/kinetics/` aliases used by the kinetics evaluator.
If both copies exist, both must match the recorded benchmark CSV identity.

## Expected format and checksums

The archived original supplementary PDF calls the additional tables Excel
files. The benchmark used CSV serializations with a descriptive first line
and a column-header second line (`skiprows=1` in the evaluator).
`upstream_references.tsv` records the SHA-256 values of those benchmark CSVs.
Current publisher downloads and Excel-to-CSV conversions remain unverified
against these identities.

For an Excel source or a different CSV serialization, obtain the matching
benchmark CSV or validate the converted table numerically. Document its
column layout, gene identities and numerical agreement before evaluation.

Original citation: Battich et al., *Sequencing metabolically labeled
transcripts in single cells reveals mRNA turnover strategies*, Science 367,
1151–1156 (2020), doi:10.1126/science.aax3072.
