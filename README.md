# RNA velocity benchmark reproduction repository

This repository curates the scripts, compact numerical results, and figure
assets needed to inspect and reproduce the RNA velocity benchmark. The
benchmark contains 20 method configurations; their sole registry is
[`methods/configurations.tsv`](methods/configurations.tsv). Large datasets,
method-result H5AD files, checkpoints, environments, logs, and exploratory
outputs remain in external storage.

The repository supports audit and reuse of most reported results, but it is not
yet a clean-environment, from-public-inputs, one-command reproduction release.

## Current status

Status as of 2026-10-09:

| Area | Repository support | Remaining boundary |
| --- | --- | --- |
| Method configurations | All 20 configurations are registered and mapped to method scripts | Reviewed environments and some upstream source revisions remain pending |
| Cell-level evaluation | Complete compact results for 14 real datasets × 20 configurations under the fixed-reference-30 contract | Large method-result inputs remain external |
| Gene-level evaluation | Reference-30 real and simulated results, repository-owned phase/comparison/summary renderers, and paper panels are available | Canonical input and method-result H5AD files remain external |
| Simulation | The three scenario configurations and a public generator are available | Runs are intentionally unseeded independent realizations, not byte-exact reconstructions of the paper inputs |
| Kinetics recovery | Frozen central-result evaluation, current Figure 6 panels/composite, and plotting entries are available; DeepKINET uses seed 123 | Rebuilding requires external H5AD files and reference tables |
| Ensemble | Mean-4 baseline and formal RowGeometricMedian-4 results are complete for all 14 datasets | Both comparators depend on the external method-result artifacts |
| Integrated summary | Current Figure 7 inputs, original editable PowerPoint and cropped author PDF are synchronized and hashed | Final PDF export uses the documented manual WPS step |
| Paper figures | Figures 1/5 and S9/S11 retain verified online PDF bytes; scientific source panels are mapped for on-demand export; Figures 4/6/7 retain validated author exports | Independent reconstruction gaps remain |

Signed-T Mean-4 and Signed-T RowGeometricMedian-4 are post-hoc comparators, not
additional registered methods. Mean-4 is the baseline and
RowGeometricMedian-4 is the formal ensemble; each is ranked independently
against the same 20 methods. See [`ensemble/README.md`](ensemble/README.md) for
the complete aggregation and ranking contract.

## Where to start

[Browse benchmark results](https://tomorin259.github.io/benchmark_repository/)
for score tables, gene phase portraits and kinetic-parameter trajectories.
Website build and source contracts are documented in [`site/README.md`](site/README.md).

| Task | Documentation and entry point |
| --- | --- |
| Inspect or run a method configuration | [`methods/README.md`](methods/README.md) |
| Reproduce Cell-level metrics and panels | [`evaluation/cell/README.md`](evaluation/cell/README.md) |
| Reproduce Gene-level metrics and phase portraits | [`evaluation/gene/README.md`](evaluation/gene/README.md) |
| Generate a new simulation realization | [`evaluation/gene/simulation/README.md`](evaluation/gene/simulation/README.md) |
| Reproduce kinetic-parameter recovery and Figure 6 panels | [`evaluation/kinetics_recovery/README.md`](evaluation/kinetics_recovery/README.md) |
| Reproduce the post-hoc ensemble | [`ensemble/README.md`](ensemble/README.md) |
| Inspect cross-experiment summaries | [`evaluation/summary/README.md`](evaluation/summary/README.md) |
| Reproduce veloVI Table S6 coordinate-sensitivity diagnostics | [`supplement/velovi/README.md`](supplement/velovi/README.md) |
| Trace manuscript panels to repository assets | [`docs/figure_map.md`](docs/figure_map.md) |
| Read shared rules, environment requirements, and historical evidence | [`docs/`](docs/) |

## Repository structure

```text
methods/                         method-specific reproduction scripts
evaluation/cell/                 Cell-level evaluation, results, and panels
evaluation/gene/                 Gene-level evaluation, simulation, and panels
evaluation/kinetics_recovery/    kinetic-parameter recovery and Figure 6 panels
evaluation/summary/              cross-experiment summaries such as Figure 7
ensemble/                        post-hoc signed-transition comparators
data/                            dataset identities and upstream access metadata
artifacts.tsv                    logical registry for large external artifacts
figures/                         paper composites and assembly assets
docs/                            shared contracts, environment requirements, figure map, and provenance
tests/                           lightweight publication-contract checks
```

Method display names are defined by `manuscript_name` in
[`methods/configurations.tsv`](methods/configurations.tsv). Dataset display
names are defined by `paper_name` in [`data/datasets.tsv`](data/datasets.tsv).
Operational IDs and file paths may retain shorter internal names.

## Reproducibility boundaries

- [`artifacts.tsv`](artifacts.tsv) records logical artifact identities and
  checksums. [`data/release_map.tsv`](data/release_map.tsv) and its verifier prepare
  publication packages; stable public URIs have not yet been published.
- [Environment documentation](docs/environments.md) records environment evidence
  and release requirements. The validated plotting lock is available; reviewed method-training
  specifications remain pending.
- Simulation generation intentionally leaves the random stream unseeded. The
  historical canonical H5AD files remain the exact paper inputs; a new run is
  an independent realization.
- [`figures/assembly/assembly_manifest.tsv`](figures/assembly/assembly_manifest.tsv)
  records available panel files, SHA256 values, and unresolved source gaps.
- CI checks lightweight publication contracts and compact results. It does not
  train all methods, install every method environment, or download every
  external artifact.

Read [`AGENTS.md`](AGENTS.md) and
[`REPOSITORY_DESIGN.md`](REPOSITORY_DESIGN.md) before changing a script or
module.

---

# RNA velocity 基准评测复现仓库

本仓库整理了检查和复现 RNA velocity 基准评测所需的脚本、紧凑数值结果和图件素材。
基准评测包含 20 个方法配置，唯一注册表是
[`methods/configurations.tsv`](methods/configurations.tsv)。大型数据集、方法结果
H5AD、模型检查点、软件环境、日志和探索性输出保存在外部存储中。

本仓库已经能够审计和复用论文中的大部分结果，但尚未达到从公开输入和干净环境
开始、一条命令完整复现全部论文结果的发布状态。

## 当前状态

截至 2026-10-09：

| 模块 | 仓库现有支撑 | 当前边界 |
| --- | --- | --- |
| 方法配置 | 20 个配置均已注册并映射到方法脚本 | 正式环境及部分上游源码版本仍待补充 |
| Cell-level 评估 | fixed-reference-30 契约下的 14 个真实数据集 × 20 个配置紧凑结果已完成 | 大型方法结果输入仍在外部存储 |
| Gene-level 评估 | reference-30 真实及模拟结果、仓库独立的相图/组图/汇总绘图入口和论文 panel 已归档 | canonical 输入和方法结果 H5AD 仍在外部存储 |
| 模拟数据 | 三个场景的配置和公开生成器已归档 | 每次运行都是不设 seed 的独立 realization，不是论文输入的 byte-exact 重建 |
| 动力学参数恢复 | central-results 评估、当前 Figure 6 单图素材、组合图和绘图入口已同步；DeepKINET 使用 seed 123 | 重建仍需外部 H5AD 和参考表 |
| 集成方法 | Mean-4 baseline 和正式 RowGeometricMedian-4 的 14 数据集结果已完成 | 两个比较器均依赖外部方法结果 artifact |
| 综合汇总 | Figure 7 当前输入表、原始可编辑 PowerPoint 和裁剪后的作者导出 PDF 已同步并记录 hash | 最终 PDF 使用文档说明的 WPS 手动导出步骤 |
| 论文图件 | Figure 1/5 与 S9/S11 保留已核对的在线 PDF 原版；科学 panel 按需导出；Figure 4/6/7 保留已验证的作者导出图 | 独立重建缺口仍有记录 |

Signed-T Mean-4 和 Signed-T RowGeometricMedian-4 是事后比较器，不属于 20 个注册
方法。Mean-4 是 baseline，RowGeometricMedian-4 是正式 ensemble；二者分别与同一组
20 个方法独立排名。完整的聚合及排名契约见
[`ensemble/README.md`](ensemble/README.md)。

## 从哪里开始

| 任务 | 文档和入口 |
| --- | --- |
| 检查或运行方法配置 | [`methods/README.md`](methods/README.md) |
| 复现 Cell-level 指标和 panel | [`evaluation/cell/README.md`](evaluation/cell/README.md) |
| 复现 Gene-level 指标和 phase portrait | [`evaluation/gene/README.md`](evaluation/gene/README.md) |
| 生成新的模拟数据 realization | [`evaluation/gene/simulation/README.md`](evaluation/gene/simulation/README.md) |
| 复现动力学参数恢复和 Figure 6 panel | [`evaluation/kinetics_recovery/README.md`](evaluation/kinetics_recovery/README.md) |
| 复现事后 ensemble | [`ensemble/README.md`](ensemble/README.md) |
| 检查跨实验汇总 | [`evaluation/summary/README.md`](evaluation/summary/README.md) |
| 复现 veloVI Table S6 坐标敏感性诊断 | [`supplement/velovi/README.md`](supplement/velovi/README.md) |
| 查询论文 panel 对应的仓库素材 | [`docs/figure_map.md`](docs/figure_map.md) |
| 阅读共用规则、环境要求和历史验证证据 | [`docs/`](docs/) |

## 仓库结构

```text
methods/                         各方法的复现脚本
evaluation/cell/                 Cell-level 评估、结果和 panel
evaluation/gene/                 Gene-level 评估、模拟数据和 panel
evaluation/kinetics_recovery/    动力学参数恢复和 Figure 6 panel
evaluation/summary/              Figure 7 等跨实验汇总
ensemble/                        事后有符号转移比较器
data/                            数据集身份和上游访问信息
artifacts.tsv                    大型外部 artifact 的逻辑注册表
figures/                         论文组合图和 assembly 素材
docs/                            跨模块契约、环境要求、figure map 和历史证据
tests/                           轻量级论文契约检查
```

方法展示名称由 [`methods/configurations.tsv`](methods/configurations.tsv) 中的
`manuscript_name` 定义；数据集展示名称由
[`data/datasets.tsv`](data/datasets.tsv) 中的 `paper_name` 定义。运行 ID 和文件路径
可以继续使用较短的内部名称。

## 复现边界

- [`artifacts.tsv`](artifacts.tsv) 已记录逻辑 artifact 身份和校验和，但稳定公开 URI
  尚未发布；[`data/release_map.tsv`](data/release_map.tsv) 与校验脚本用于准备发布包。
- [环境说明](docs/environments.md) 记录环境证据与发布要求，已提供验证过的绘图 lock；方法训练环境
  specification 仍待补充。
- 模拟数据生成器刻意不设置随机种子。历史 canonical H5AD 仍是论文使用的精确输入；
  新运行得到的是独立随机 realization。
- [`figures/assembly/assembly_manifest.tsv`](figures/assembly/assembly_manifest.tsv)
  记录现有 panel 文件、SHA256 和未解决的 source 缺口。
- CI 只检查轻量级论文契约和紧凑结果，不会训练全部方法、安装所有方法环境或下载
  全部外部 artifact。

修改脚本或模块前，请先阅读 [`AGENTS.md`](AGENTS.md) 和
[`REPOSITORY_DESIGN.md`](REPOSITORY_DESIGN.md)。
