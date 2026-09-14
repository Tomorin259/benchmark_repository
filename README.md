# RNA velocity benchmark reproduction repository

This repository curates the code and compact numerical outputs needed to
understand and reproduce the RNA velocity benchmark. The benchmark has 20
method configurations, registered only in
methods/configurations.tsv. Large datasets, method result H5AD files,
checkpoints, environments, logs, and exploratory outputs remain in external
storage.

The two post-hoc comparators are separate from those 20 methods: Signed-T
Mean-4 is a baseline and Signed-T RowGeometricMedian-4 is the formal ensemble.
Each is ranked against
the 20 methods separately; the two comparators are not ranked against each
other. Mean-4 is the arithmetic mean of the four signed transition values on
each common edge. RowGeometricMedian-4 takes the Euclidean geometric median of
the four complete signed off-diagonal transition rows independently for each
source cell, without post-hoc L1 renormalization. For each dataset and metric,
each variant is ranked independently
against the 20 Cell methods: rank is 1 plus the number of strictly higher method
scores, with exact method ties reported using competition-rank tie counts; the
ensemble is participant 21 and variants are never candidates for one another.
Their implementations are in ensemble/, using the portable fixed-neighbor
assets in evaluation/cell/reference_neighbors/. Those assets install the
audited candidate graph and indices, not their build-time PCA coordinates;
cosines are still computed in each input method's native Ms/velocity space.
Formal fixed-neighbor all-14 ensemble results are
complete and registered in `artifacts.tsv`. The canonical scores, summaries, and
independent rankings are [`ensemble/results/scores.csv`](ensemble/results/scores.csv),
[`ensemble/results/summary.md`](ensemble/results/summary.md),
[`ensemble/results/rankings_vs_20.csv`](ensemble/results/rankings_vs_20.csv), and
[`ensemble/results/ranking_summary.md`](ensemble/results/ranking_summary.md).

## Current status

As of 2026-09-14, the repository has a complete compact Cell-level result
package for the 14 real datasets and 20 method configurations, including the
fixed-reference-30 contract and the formal Mean-4 and RowGeometricMedian-4
all-14 ensemble results. Gene-level and kinetics-recovery outputs are curated
result packages with analysis and plotting entry points, but they are not yet
an end-to-end rebuild from public raw inputs. The integrated summary contains
the compact Figure 6 input tables, while its final renderer is still missing.

The repository is therefore suitable for result inspection and contract-level
auditing, not yet for a clean-environment, from-scratch reproduction release.
In particular, `envs/` contains release requirements but no reviewed Conda
specification or lock file, while `data/` and `artifacts.tsv` record dataset
and artifact identities but the external URIs are not published. The three
simulation scenarios now have an in-repository generator and audited scenario
configuration. The generator intentionally leaves the random stream unseeded,
so each run is an independent realization; the historical canonical H5AD
files remain the exact paper inputs and are not claimed to be byte-exact
regenerations. The figures assembly manifest records these limits and the
remaining missing panel sources explicitly. The presence of a directory must
not be read as full all-20 benchmark completion.

## Repository path

```text
methods/                         method-specific reproduction scripts
envs/                            environment release requirements; specifications pending
evaluation/cell/                 cell-level evaluation and panels
evaluation/gene/                 gene-level evaluation and simulation
evaluation/kinetics_recovery/  kinetic-parameter recovery and its plotting module
evaluation/summary/              cross-experiment summaries such as Figure 6
ensemble/                        optional post-hoc signed-transition utility
data/                            dataset metadata and access information
artifacts.tsv                    portable index of external artifacts
figures/                         confirmed paper-level composite PDFs
docs/                            cross-module documentation and figure map
tests/                           lightweight publication-contract checks
```

There is intentionally no shared top-level preprocessing package or results
directory. The top-level figures directory contains only confirmed paper-level
composite PDFs; panel assets stay in the evaluation module that owns the
scientific question.

Method names shown in public tables use veloVI spelling. Dataset IDs identify
the biological input and evaluation panel, not methods or ensemble variants.

## Running kinetics-recovery

The kinetics-recovery module consumes saved method-result H5AD files and the
scEU-seq reference tables. It does not train methods or prepare the benchmark
datasets.

```bash
python evaluation/kinetics_recovery/metrics.py \
  --artifact-root /path/to/artifacts \
  --check-only

python evaluation/kinetics_recovery/metrics.py \
  --artifact-root /path/to/artifacts \
  --output /path/to/metrics

python evaluation/kinetics_recovery/plotting/figure5.py \
  --metrics /path/to/metrics \
  --output /path/to/figures
```

The evaluator extracts each method's native parameter keys at its input
boundary, applies the fixed RPE1/organoid processing rules, calculates
per-gene and summary metrics, and writes compact CSV files. Plotting consumes
only that metrics bundle.

## Counts-specific inputs

Only methods that genuinely require raw counts select the corresponding
`_counts.h5ad` artifact and perform its method-specific checks. DeepKINET and
cell2fate document the four special datasets, cell/gene alignment, finite and
non-negative count checks, integer-valued checks, and lossless dtype conversion
inside their own method scripts. DeepKINET additionally records the seed-123
deterministic benchmark setting used by the formal result.
Methods that do not require counts do not use this path.

## Public artifact index

The public-facing artifact index is the root artifacts.tsv; it must not expose server-specific working paths.

Read AGENTS.md and REPOSITORY_DESIGN.md before changing a script or module.

---

# RNA velocity 基准评测复现仓库

本仓库整理了理解和复核 RNA velocity 基准评测所需的代码与紧凑数值结果。
基准评测包含 20 个方法配置，并统一登记在
[`methods/configurations.tsv`](methods/configurations.tsv) 中。大型数据集、
方法结果 H5AD、模型检查点、软件环境、日志和探索性输出保存在外部存储中。

两个事后比较器不属于上述 20 个方法：Signed-T Mean-4 是基线，Signed-T
RowGeometricMedian-4 是正式集成方法。二者分别与 20 个方法独立排名，彼此不参与
对方的排名。Mean-4 对公共边上的四个有符号转移值取算术平均；
RowGeometricMedian-4 则针对每个源细胞，分别计算四个完整有符号非对角转移行的
欧氏几何中位数，并且不进行事后 L1 重新归一化。对每个数据集和指标，每个比较器
都独立与 20 个 Cell 方法排名：名次等于 1 加上严格高于该分数的方法数；完全相同
的方法分数采用竞赛排名并记录并列数量。集成方法作为第 21 个参与者，两个比较器
不会互相作为候选方法。

实现位于 [`ensemble/`](ensemble/)，并使用
[`evaluation/cell/reference_neighbors/`](evaluation/cell/reference_neighbors/)
中的可移植 fixed-neighbor 资源。这些资源保存经过审计的候选图和索引，而不保存
构建时的 PCA 坐标；余弦值仍在每个输入方法自身的 `Ms`/velocity 空间中计算。
14 个数据集的正式 fixed-neighbor 集成结果已经完成并登记在 `artifacts.tsv`。
正式分数、汇总和独立排名分别位于：

- [`ensemble/results/scores.csv`](ensemble/results/scores.csv)
- [`ensemble/results/summary.md`](ensemble/results/summary.md)
- [`ensemble/results/rankings_vs_20.csv`](ensemble/results/rankings_vs_20.csv)
- [`ensemble/results/ranking_summary.md`](ensemble/results/ranking_summary.md)

## 当前状态

截至 2026-09-14，本仓库已经形成完整、紧凑的 Cell-level 结果包，覆盖 14 个
真实数据集和 20 个方法配置，包括 fixed-reference-30 契约，以及 Mean-4 和
RowGeometricMedian-4 在全部 14 个数据集上的正式集成结果。Gene-level 和
kinetics-recovery 已整理为包含分析及绘图入口的结果包，但尚不能从公开原始输入
端到端重建。Integrated summary 已保存 Figure 6 的紧凑输入表，但正式 renderer
仍然缺失。

因此，本仓库目前适合结果检查和契约级审计，但还不是可以在干净环境中从头完整
复现的发布版本。具体而言，`envs/` 已记录发布要求，但尚无经过审核的 Conda
环境文件或 lock file；`data/` 和 `artifacts.tsv` 已登记数据集及 artifact 身份，
但外部 URI 尚未发布。三个模拟场景现在已经具有仓库内生成器和经过审计的场景
配置。生成器刻意不设置随机种子，因此每次运行都是独立随机 realization；历史
canonical H5AD 仍作为论文使用的精确输入，本仓库不声称新生成结果可以与其
byte-exact 一致。图件 assembly manifest 明确记录了这些边界和仍缺失的 panel
来源。目录存在本身不代表该模块已经完成全部 20 个方法的端到端复现。

## 仓库结构

```text
methods/                         各方法的复现脚本
envs/                            环境发布要求；正式环境文件仍待补充
evaluation/cell/                 细胞级评估与图件
evaluation/gene/                 基因级评估与模拟数据
evaluation/kinetics_recovery/    动力学参数恢复及绘图模块
evaluation/summary/              Figure 6 等跨实验汇总
ensemble/                        事后有符号转移集成工具
data/                            数据集元数据与访问说明
artifacts.tsv                    外部 artifact 的可移植索引
figures/                         已确认的论文级组合 PDF
docs/                            跨模块文档和论文图映射
tests/                           轻量级论文契约检查
```

仓库刻意不设置共享的顶层预处理包或统一 `results/` 目录。顶层 `figures/` 只保存
已确认的论文级组合 PDF；panel 素材保存在负责相应科学问题的 evaluation 模块中。

公开表格中的方法名称统一使用 `veloVI` 拼写。数据集 ID 表示生物学输入和评估
panel，而不是方法或集成变体。

## 运行 kinetics-recovery

kinetics-recovery 模块读取已经保存的方法结果 H5AD 和 scEU-seq 参考表；它不负责
训练方法或准备基准数据集。

```bash
python evaluation/kinetics_recovery/metrics.py \
  --artifact-root /path/to/artifacts \
  --check-only

python evaluation/kinetics_recovery/metrics.py \
  --artifact-root /path/to/artifacts \
  --output /path/to/metrics

python evaluation/kinetics_recovery/plotting/figure5.py \
  --metrics /path/to/metrics \
  --output /path/to/figures
```

评估器在输入边界提取每种方法原生的参数 key，应用固定的 RPE1/organoid 处理
规则，计算逐基因指标和汇总指标，并写出紧凑 CSV。绘图模块只读取该 metrics
结果包。

## 需要原始 counts 的输入

只有确实需要原始 counts 的方法才会选择对应的 `_counts.h5ad` artifact，并执行
方法自身的输入检查。DeepKINET 和 cell2fate 在各自脚本中记录了四个特殊数据集、
细胞/基因对齐、有限且非负的 count 检查、整数值检查和无损 dtype 转换。
DeepKINET 还记录了正式结果所使用的 seed-123 确定性基准设置。不需要 counts 的
方法不会使用这一路径。

## 公开 artifact 索引

面向公开发布的 artifact 索引是仓库根目录的
[`artifacts.tsv`](artifacts.tsv)，其中不得包含服务器专用工作路径。

修改脚本或模块前，请先阅读 [`AGENTS.md`](AGENTS.md) 和
[`REPOSITORY_DESIGN.md`](REPOSITORY_DESIGN.md)。
