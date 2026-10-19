# Third-party notices

The root [MIT License](LICENSE) covers original benchmark source code.
The retained upstream files, adapted source portions and patches listed below
retain their applicable upstream terms. Installed method packages retain their
own licenses; calling their APIs does not relicense those packages.

## MIT-licensed components

### GRAVITY BATC

- Local: [Cell BATC](evaluation/cell/batc.py) and the BATC portions of
  [Ensemble metrics](ensemble/metrics.py).
- Source: [GRAVITY BATC](https://github.com/CSUBioGroup/GRAVITY/blob/6e0f2b5af8310fa469f8744df9a2e4ffbecbeeb0/gravity/analysis/batc.py).
  This revision is the source-comparison reference, not a recovered historical
  training or execution pin.
- Changes: benchmark interfaces, handling of NaN branch scores and diagnostic
  outputs. BATC's upstream attribution is retained.

### TSvelo patch

- Local: [TSvelo source patch](methods/source_patches/tsvelo.patch).
- Source: [TSvelo](https://github.com/lijc0804/TSvelo/tree/6e0ab8264e02aa511576278219fd7c1c65efec80).
- Changes: benchmark input, resource and output paths and runtime options.
  The upstream copyright line, including its brackets, is retained below.

### DeepKINET patch

- Local: [DeepKINET source patch](methods/source_patches/deepkinet.patch).
- Source: [DeepKINET](https://github.com/3254c/DeepKINET/tree/fd9cb5de02e9139322fdef12d1629cf9fc2929e5).
- Changes: package-relative imports and checkpoint-path configuration.

The following original copyright notices and MIT terms apply to the respective
upstream components above:

```text
Copyright (c) 2025 CSUBioGroup
Copyright (c) [2026] [Jiachen Li]
Copyright (c) 2023 Chikara Mizukoshi

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

## TIvelo 0.1.4 — BSD 2-Clause

- Local: `tivelo_frozen_scores` in
  [Cell metric preparation](evaluation/cell/prepare_cell_metrics.py).
- Source: the CBDir-HD and transition-score calculations in
  [TIvelo metrics](https://github.com/cuhklinlab/TIVelo/blob/v0.1.4/tivelo/utils/metrics.py).
- Changes: combined cross-environment QA calculation and benchmark result
  aggregation. The official scoring worker separately calls the installed package.
- The distributed LICENSE contains the BSD 2-Clause terms below, although the
  installed package metadata labels its license as MIT; that discrepancy does
  not replace the distributed license notice.

```text
BSD 2-Clause License

Copyright (c) 2024, MIAO Jishuai

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

## MURK screening sources — GPL version 3

- Source: [Barile et al. screening notebook](https://github.com/mebarile/Gata1_Erythroid_kinetics/blob/1e0de6338a85a559478627927c24bdc72d2d8919/Figure_3/search_murk_mouse_fig_3d_and_rank.ipynb).
- Retained original:
  [MURK notebook](evaluation/gene/selection/upstream/search_murk_mouse_fig_3d_and_rank.ipynb).
- Historical commented source excerpt:
  [Gastrulation excerpt](evaluation/gene/selection/historical/scripts/gastrulation_murk_commented_excerpt.py).
- The original notebook remains unchanged. The historical excerpt records
  upstream code as comments with benchmark context. Source identities and
  historical snapshots are documented in
  [selection provenance](evaluation/gene/selection/provenance.json).

The upstream notebook, retained excerpt and any upstream-derived portions remain
subject to the upstream GPL version 3 terms, whose full text is retained in
[LICENCE.txt](evaluation/gene/selection/upstream/LICENCE.txt). The root MIT grant
does not replace those terms. Independently written branching-selection code
has its method lineage described in the [selection guide](evaluation/gene/selection/README.md).

## TFvelo patch — CC BY 4.0

- Local: [TFvelo source patch](methods/source_patches/tfvelo.patch).
- Source: [TFvelo software archive, version 1.0](https://zenodo.org/records/10499309),
  whose recorded creator is `lijc0804`, and the
  [pinned source](https://github.com/xiaoyeye/TFvelo/tree/ec6cdb940af94f02fe32c8fcdb98494cdb4beb96).
- The official archive identifies its license as
  [Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/).
  The source files touched by this patch match that release.
- Changes: input, resource and output paths and analysis command-line options.
  Source identities and patch hashes are recorded in
  [the method source manifest](methods/source_manifest.json).

Attribution is retained to the TFvelo project and official software archive.
External ENCODE/ChEA resources retain their respective source terms and are not
relicensed by this code notice.

## MultiVelo 0.1.0 diagnostic control — BSD 3-Clause

- Local: [candidate-gate patch](supplement/multivelo/results/candidate_gate_only.patch)
  and the in-memory modification performed by
  [the diagnostic control](supplement/multivelo/run_multivelo_controls.py).
- Source: `dynamical_chrom_func.py` in
  [MultiVelo](https://github.com/welch-lab/MultiVelo), package version 0.1.0.
- Change: retain repression candidates by disabling one candidate-exclusion
  condition for the diagnostic comparison.

```text
BSD 3-Clause License

Copyright (c) 2021, Welch Lab
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

## cellDancer simulator and sampling

The three simulated benchmark scenarios were generated using the
`celldancer.simulation.simulate` interface in cellDancer 1.1.7; see the
[simulation workflow](evaluation/gene/simulation/README.md).

The [phase-portrait API](evaluation/gene/plotting/phase_portrait.py) was written
for this benchmark. Its representative-cell sampling strategy was inspired by
cellDancer 1.1.7. The existing
[sampling-source notice](evaluation/gene/plotting/THIRD_PARTY_NOTICES.md)
retains the upstream BSD 3-Clause notice. This acknowledgment does not attribute
the complete plotting API to cellDancer.

## DeepVelo_scv

DeepVelo_scv follows the official DeepVelo Figure2 workflow.
