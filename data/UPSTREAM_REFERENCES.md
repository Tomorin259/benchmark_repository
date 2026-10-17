# Obtain the original kinetic-reference tables

This distribution provides 25 mirrored scientific files and records three
additional upstream-only reference dependencies. It does not host Battich's
original Tables S1/S3/S4 or reproduce their reference trajectories.

1. Visit https://doi.org/10.1126/science.aax3072 and open Supplementary Materials.
2. Obtain Tables S1, S3 and S4 through the original publisher/author access route.
3. Place the benchmark-compatible CSV files at the paths in
   `upstream_references.tsv`, using the layout required by your evaluator.
4. From this repository root, run `python data/check_references.py --artifact-root /path/to/artifacts`.
   In the HF version folder, run `python check_references.py --artifact-root /path/to/artifacts` instead.

There are two documented layouts: velocity-root paths used by existing
launchers and `references/kinetics/` aliases used by the kinetics evaluator.
If both copies exist, both must match the recorded benchmark CSV identity.

## Format and verification boundary

The archived original supplementary PDF calls the additional tables Excel
files. The benchmark used CSV serializations with a descriptive first line
and a column-header second line (`skiprows=1` in the evaluator). This release
records hashes of those exact benchmark CSVs; they are not claimed to be
hashes of current publisher downloads. Automated access to the publisher
was blocked during preparation, so no guessed direct file URLs or unverified
Excel-to-CSV conversion are supplied.

If the obtained file is Excel, or its CSV serialization differs, do not merely
rename it or suppress a failed check. Obtain the matching CSV or establish
and document a conversion with numerical validation before reporting an
exact reproduction. This release does not claim one-command retrieval or
byte-identical regeneration of a differently formatted upstream file.

Original citation: Battich et al., *Sequencing metabolically labeled
transcripts in single cells reveals mRNA turnover strategies*, Science 367,
1151–1156 (2020), doi:10.1126/science.aax3072.
