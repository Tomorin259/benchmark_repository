# RNA velocity benchmark 项目导航与图件溯源

> 最后核验：2026-08-31
> 服务器项目根目录：`/mnt/second19T/zhaozelin/velocity`
> 论文主文件：Overleaf 项目中的 `merged.tex`

## 1. 本文用途

本文是本项目的快速入口，供后续研究者或 agent 在开始工作前建立统一上下文。它主要回答四个问题：

1. 项目各目录负责什么；
2. 论文主图和补充图实际由哪些结果、脚本和中间文件生成；
3. 哪些目录属于当前结果，哪些属于历史版本或诊断实验；
4. 修改、重跑或清理文件之前需要先检查什么。

本文记录的是已经从 Overleaf 最终 PDF 反向核验过的文件关系。目录名只能作为线索，不能替代最终图件、数值、manifest 或文件哈希的核验。

## 2. 项目定位

本项目以 RNA velocity benchmark 为主线，对 20 个方法进行细胞层面和基因层面的系统评价，主要包括：

- 14 个真实数据集上的细胞轨迹评价；
- 模拟数据上的基因级方向与相似性评价；
- 真实数据中 MURK genes 和 branching genes 的动力学表现；
- RPE1 和 organoid scEU-seq 数据上的动力学参数恢复评价；
- benchmark 之后附带提供的 signed-transition ensemble。

`transition_ensemble` 不是第 21 个 RNA velocity 模型，而是一个不学习方法权重、只面向 cell-level signed transition operator 的事后集成工具。

## 3. 顶层目录

| 目录 | 主要职责 | 使用注意事项 |
|---|---|---|
| `data/` | 原始或预处理数据 | 不应根据绘图需求直接改写 |
| `methods/` | 各 RNA velocity 方法的代码、环境和原始运行输出 | `methods/*/results*` 不一定是论文成图使用的 canonical result |
| `results/` | 经整理后的方法结果，很多细胞层面成图从这里读取 | 具体文件必须以绘图脚本或 manifest 为准 |
| `indicators/cell/` | 细胞层面指标计算及汇总结果 | 包括 CBDir、BATC、Consistency、ICCoh 和 TransCosine 等 |
| `indicators/gene/` | 模拟、真实 branching/MURK 及 kinetic-parameter 分析 | 包含若干历史复现和诊断目录，不能只按文件名判断版本 |
| `plot/` | 论文图的绘图脚本、当前输出和历史输出 | 最终 PDF 多数经过人工拼版，需按 panel 追溯 |
| `MURK/` | MURK gene 相关候选集合与辅助材料 | 主要作为 gene-level 分析的上游依据 |
| `transition_ensemble/` | signed-transition ensemble 的独立 Git 仓库 | 修改前必须检查分支、HEAD 和工作区状态 |
| `archive/` | 项目级归档 | 不默认参与当前成图 |
| `plot/archive_old_versions_20260613/` | 2026-06-13 前后的绘图与结果备份 | **不能整体删除：Fig. 5A 仍使用其中五张流线图** |

关键原则：

> 最终图件的直接输入通常位于 `results/`、`indicators/` 或 `plot/`；`methods/` 中的脚本和结果通常是更上游的来源。

## 4. 论文与最终图件

论文以 Overleaf 中的 `merged.tex` 为唯一主线。主要图件文件为：

- `main_figs/figOne.pdf` 至 `main_figs/figSix.pdf`；
- `supplement/suppl_figs/suppFigOne.pdf` 至 `suppFigTwelve.pdf`。

最终 PDF 往往是从多个 PNG、PDF 或 SVG 手工组合而成。以下记录的是已经确认的实验和素材来源，而不一定是一个能够独立生成完整最终 PDF 的脚本。

## 5. 主图溯源

### Figure 1

Figure 1 是概念性 benchmark 流程图。

- 未在 `velocity/` 下找到完整的自动生成脚本或可编辑源文件；
- Methods、数据集表和方法配置是其内容依据，不是可执行的生成链路；
- 当前应将其视为外部手工制作并上传至 Overleaf 的图件。

因此，在没有找到原始 PPT、AI、SVG 或其他编辑源之前，不应声称 Figure 1 可以从仓库自动复现。

### Figure 2：14 个真实数据集的 cell-level benchmark

#### Panel A–B：流线图与参考方向

主要脚本：

```text
plot/unified_velocity_stream/plot_unified_streams.py
plot/unified_velocity_stream/build_velocity_montages_3x7.py
```

当前最终素材：

```text
plot/unified_velocity_stream/montages/gastrulation_20methods_3x7_BATC.png
plot/unified_velocity_stream/montages/DentateGyrus_lamanno_20methods_3x7_CBDir.png
```

精确结果路径和细胞对齐情况记录在：

```text
plot/unified_velocity_stream/results/render_manifest.csv
```

该 manifest 是判断 20 个方法实际读取了哪个 H5AD 的首要依据。不要根据 `methods/*` 下的同名结果猜测输入文件。

#### Panel C：五个 cell-level 指标

绘图脚本：

```text
evaluation/cell/plots.py
```

直接输入：

```text
evaluation/cell/results/all_metrics.csv
evaluation/cell/results/transprob.csv
```

最终组合图：

```text
evaluation/cell/figures/final_seed123/fig2_cell_metrics_5panel.pdf
evaluation/cell/figures/final_seed123/fig2_cell_metrics_5panel.png
evaluation/cell/figures/final_seed123/fig2_cell_metrics_5panel.svg
```

相关计算入口包括：

```text
evaluation/cell/metrics.py
evaluation/cell/transprob.py
```

运行时结果映射位于 `evaluation/cell/result_inputs.tsv`。它与
`manifests/cell_result_inputs.tsv` 的流线图审计映射分开维护；例如
DeepVelo 的 cell metrics 使用 metric-result 文件族，而统一流线图可以使用
`results_with_vu` 文件族。

明确未使用的旧版：

```text
plot/archive_old_versions_20260613/before_cell_pdf_rerun_20260613_034316/
plot/cell_metrics_5panel/plot_cell_metric_panels.py.bak_leftmargin_20260702_1348
```

### Figure 3：模拟数据的 gene-level benchmark

#### Panel A：模拟任务示意

```text
plot/to_use_sim/overall/plot_new_genes.py
```

#### Panel B–C：方法级 phase portrait/velocity 对比

```text
plot/to_use_sim/details/plot_new_method_comparisons.py
```

直接读取整理后的模拟结果：

```text
indicators/gene/simulation/results_for_analyze/
```

`methods/*/*simulation.py` 属于上游方法运行过程，不是 Figure 3 的直接绘图输入。

#### Panel D：CBDir 与 CosSim 汇总

```text
plot/to_use_sim/details/plot_new_score_summary.py
```

直接读取三个当前 score 文件：

```text
indicators/gene/simulation/branching_analyze/results_new/score_rerun_results.txt
indicators/gene/simulation/murk_analyze/results_new/score_rerun_results.txt
indicators/gene/simulation/sim_all_analyze/results_new/score_rerun_results.txt
```

当前图输出：

```text
plot/to_use_sim/details/results_new/
```

明确的旧版输出：

```text
plot/to_use_sim/details/results_new_before_deepvelo_sim_rerun_20260613/
indicators/gene/simulation/branching_analyze/results_new_before_deepvelo_sim_rerun_20260613/
indicators/gene/simulation/murk_analyze/results_new_before_deepvelo_sim_rerun_20260613/
indicators/gene/simulation/sim_all_analyze/results_new_before_deepvelo_sim_rerun_20260613/
```

DeepVelo 的当前模拟输入为：

```text
indicators/gene/simulation/results_for_analyze/deepvelo/deepvelo_{all,branching,murk}.h5ad
```

同目录下带 `20260613_before_simulation_rerun` 的文件是重跑前备份。

2026-07-14 生成的下列 TIvelo 文件是 `filter_genes` 诊断实验，并未替换 Figure 3 的 canonical TIvelo 输入：

```text
indicators/gene/simulation/results_for_analyze/tivelo_filter_genes_*/
indicators/gene/simulation/results_for_analyze/tivelo_filter_genes_pair_seed20260714_comparison/
```

### Figure 4：真实数据中的 branching/MURK gene 分析

主要绘图与计算脚本：

```text
plot/touse_murk/murk_analyze.py
plot/gene_to_use/branching_total_results_from_final_genes.py
```

当前输出：

```text
plot/final_genes/
plot/gene_to_use/branching_combined_from_final_genes.png
```

`murk_analyze.py` 对每个数据集读取九个 gene-level-compatible 方法结果：

1. cellDancer；
2. DeepVelo；
3. DeepVelo 2022/scVelo reconstruction；
4. GraphVelo-scVelo；
5. MultiVelo；
6. scVelo dynamical；
7. TIvelo；
8. UniTVelo with `velocity_u`；
9. veloVI with `velocity_u`。

当前 Figure 4 中可读出的 Gastrulation/Abcg2 分数与 `plot/final_genes/gastrulation/score.txt` 一致，不与旧版一致。

明确的旧版：

```text
plot/archive_old_versions_20260613/before_gene_cbdir_rerun_20260613_013621/final_genes/
```

版本注意：当前 `murk_analyze.py` 的修改时间晚于 `plot/final_genes/` 的生成时间。由于 `velocity/` 根目录没有 Git 历史，目前只能确认它是当前维护脚本，不能仅凭时间戳证明它与已上传图件所用代码逐字节一致。

### Figure 5：scEU-seq kinetic-parameter benchmark

Figure 5 必须按 panel 追溯，不能只给一个目录。

#### Panel A：RPE1 流线图

最终使用以下五张素材：

```text
plot/archive_old_versions_20260613/rpe1/scvelo__stream2_dynamo.png
plot/archive_old_versions_20260613/rpe1/scvelo__stream2_celldancer.png
plot/archive_old_versions_20260613/rpe1/scvelo__stream2_deepvelo.png
plot/archive_old_versions_20260613/rpe1/scvelo__stream2_deepkinet.png
plot/archive_old_versions_20260613/rpe1/scvelo__stream2_tsvelo.png
```

这是 Figure 5 中唯一直接使用 `archive_old_versions_20260613` 的部分。它表示最终选中的静态素材被移动进了 archive，并不表示整个 Figure 5 使用旧实验结果。

#### Panel B：RPE1 alpha/gamma trajectory curves

```text
evaluation/kinetics_recovery/figures/final_seed123/rpe1_alpha_trajectories.png
evaluation/kinetics_recovery/figures/final_seed123/rpe1_gamma_trajectories.png
```

生成脚本：

```text
evaluation/kinetics_recovery/plots.py
```

#### Panel C：RPE1 alpha/gamma heatmap

```text
evaluation/kinetics_recovery/figures/final_seed123/rpe1_alpha_gamma_heatmap.png
```

它与下列文件的 SHA-256 完全一致：

```text
evaluation/kinetics_recovery/figures/final_seed123/rpe1_alpha_gamma_heatmap.png
```

#### Panel D：kinetic-class Venn diagrams

- 左侧 cellDancer 三图来自 `indicators/gene/scEU_rpe1_analyze/venn_cellDancer.png`；
- 右侧 DeepVelo 三图来自 `figure5_candidates/rpe1_deepvelo_venn.png`；
- DeepVelo 的计数与 `figure5_candidates/manifest.json` 一致。

#### Panel E：organoid 流线图

生成脚本：

```text
plot/figure5_panel_e_original_style/plot_organoid_original_celldancer_style.py
```

当前素材：

```text
evaluation/kinetics_recovery/figures/final_seed123/organoid_streams/celldancer_organoid_stream_title_removed.png
evaluation/kinetics_recovery/figures/final_seed123/organoid_streams/DeepKINET_organoid_stream_title_removed.png
evaluation/kinetics_recovery/figures/final_seed123/organoid_streams/deepvelo_organoid_stream_title_removed.png
```

DeepKINET 的最终 H5AD 不预存 velocity graph；重画时在流线图前按
`DeepKINET_velocity` 显式重算 graph。

#### Panel F：organoid branch-specific boxplots

```text
evaluation/kinetics_recovery/figures/final_seed123/organoid_branch1_boxplot.png
evaluation/kinetics_recovery/figures/final_seed123/organoid_branch2_boxplot.png
```

上游 DeepVelo kinetic 数据由以下统一入口生成：

```text
evaluation/kinetics_recovery/metrics.py
```

### Case 3 kinetic-parameter evaluation

Case 3 直接读取既有方法结果，不再复制数据准备或方法训练步骤。入口为：

~~~text
evaluation/kinetics_recovery/config.json
evaluation/kinetics_recovery/metrics.py
~~~

case3.json 只维护参考表、结果 H5AD 路径、参数 layer 映射和必要的
方法特异提取规则；
metrics.py 负责 RPE1 与 organoid 的轨迹处理、逐基因 Spearman/cosine、
component means 和 Figure 6 的 equal-system summary。cellDancer、
DeepVelo 与 DeepKINET 直接引用现有 H5AD；Dynamo 和 TSvelo 直接引用
各自专用 H5AD：Dynamo 从带标记结果的 sparse cell_wise_alpha 提取，
TSvelo 从两个 lineage H5AD 的 alpha 合并提取。二者的历史 CSV 仅用于
数值回归，不是运行输入。plots.py 只读取 metrics 输出并生成
Figure 5/6 的定量 review panels 与绘图表，不重复读取或处理上游数据。

### Figure 6：方法属性与跨层级表现汇总

当前直接数值来源：

```text
evaluation/summary/results/figure6_summary/fig2_gene_metric_wps_values.csv
evaluation/summary/results/figure6_summary/fig2_gene_metric_summary.csv
evaluation/summary/results/figure6_summary/fig2_gene_metric_dataset_means.csv
evaluation/summary/results/figure6_summary/fig2_kinetic_component_means.csv
```

当前可找到的重建脚本：

```text
evaluation/summary/README.md
```

该脚本负责 kinetic Spearman 的 alpha/gamma 分离、跨 system 汇总和展示值 `+1` 转换，但它不是整张 Figure 6 的完整生成程序。

已知限制：

- 最终 Figure 6 应为人工组合；
- 没有找到从全部底层指标生成 `fig2_kinetic_component_means.csv` 的脚本；
- 脚本和 CSV 中保留的 `fig2` 名称是内部历史命名，论文中对应 Figure 6；
- 当前未发现这些四个 CSV 的同名历史副本。

## 6. Supplementary Figures

### Figures S1–S8

Figures S1–S8 与 Figure 4 共用当前 `plot/final_genes/` 流程：

| 补充图 | 数据集目录 |
|---|---|
| S1 | `plot/final_genes/DentateGyrus_lamanno/` |
| S2 | `plot/final_genes/retina/` |
| S3–S4 | `plot/final_genes/gastrulation/` |
| S5 | `plot/final_genes/capo_pre/` |
| S6 | `plot/final_genes/bonemarrow/` |
| S7 | `plot/final_genes/mousebrain/` |
| S8 | `plot/final_genes/pancreas/` |

S7 中 TIvelo 对 Adcy1、Cacna1e、Flrt2 和 Gria2 的显示值约为 `0.46/0.38/0.41/0.50`，逐项对应当前 `mousebrain/score.txt`，不对应归档版本。这与 Figure 4 的数值核验共同表明，最终 Figure 4 和 S1–S8 使用的是 2026-06-13 指标重跑后的 `final_genes` 系列。

### Figures S9–S12：signed-transition ensemble

核心代码：

```text
transition_ensemble/ensemble.py
transition_ensemble/geometric_median.py
transition_ensemble/benchmark.py
transition_ensemble/plot_results.py
```

当前 14 数据集结果：

```text
transition_ensemble/results/row_geometric_median_all14/scores.csv
transition_ensemble/results/row_geometric_median_all14/mean4/
transition_ensemble/results/row_geometric_median_all14/row_geometric_median4/
```

论文补充图的当前含义为：

| 补充图 | 内容 | 直接来源 |
|---|---|---|
| S9 | signed-transition ensemble 工作流示意图 | 论文侧 `scripts/build_ensemble_supplement.py` |
| S10 | Mean-4 与 RowGeometricMedian-4 的 14 数据集配对比较 | `scores.csv` |
| S11 | GeoMedian-4 与原 20 方法结果的背景比较 | ensemble `figures/` 中两张背景箱线图 |
| S12 | DentateGyrus 与 organoid 的 Mean-4/GeoMedian-4 代表性流线图 | `mean4/` 与 `row_geometric_median4/` 中四张 stream PNG |

S11 的当前直接素材为：

```text
transition_ensemble/results/row_geometric_median_all14/figures/Benchmark_Fig2_CBDir_with_row_geometric_median.png
transition_ensemble/results/row_geometric_median_all14/figures/Benchmark_Fig2_BATC_with_row_geometric_median.png
```

同目录中带 `with_GeoMedian4` 的 PDF/PNG/SVG 也是当前结果素材，但 2026-08-31 核验的论文构建脚本仍直接读取上述两个历史内部命名的 PNG；发布时应统一命名并同步修改脚本。

S12 的四张直接素材为：

```text
transition_ensemble/results/row_geometric_median_all14/mean4/dentategyrus.stream.png
transition_ensemble/results/row_geometric_median_all14/row_geometric_median4/dentategyrus.stream.png
transition_ensemble/results/row_geometric_median_all14/mean4/organoid.stream.png
transition_ensemble/results/row_geometric_median_all14/row_geometric_median4/organoid.stream.png
```

当前 S9–S12 的论文侧构建脚本为：

```text
/Users/tomorin/Downloads/test1_1/scripts/build_ensemble_supplement.py
```

该脚本直接生成：

```text
/Users/tomorin/Downloads/test1_1/supplement/suppl_figs/suppFigNine.pdf
/Users/tomorin/Downloads/test1_1/supplement/suppl_figs/suppFigTen.pdf
/Users/tomorin/Downloads/test1_1/supplement/suppl_figs/suppFigEleven.pdf
/Users/tomorin/Downloads/test1_1/supplement/suppl_figs/suppFigTwelve.pdf
```

截至 2026-08-31，远端 `transition_ensemble` 位于 `main` 分支，工作区干净，当前 HEAD 为：

```text
de47402ea8c5625f765c595b279d8d53a4594cec
```

14 数据集的 28 个 CBDir/BATC 分数已经与当前论文表核对，最大绝对差为 `1.11e-16`。引用或重跑 S9–S12 前仍应执行：

```bash
cd /mnt/second19T/zhaozelin/velocity/transition_ensemble
git status --short
git log -5 --oneline
```

## 7. 当前结果与历史结果的判断规则

### Cell-level

1. 最终分数先看 `indicators/cell/all_metrics/*.csv`；
2. 流线图所用 H5AD 先看 `plot/unified_velocity_stream/results/render_manifest.csv`；
3. 不使用 `before_cell_pdf_rerun_*` 中的旧 PDF；
4. 不根据方法目录中最接近的文件名推断 canonical result。

### Simulation gene-level

1. 绘图输入以 `indicators/gene/simulation/results_for_analyze/` 为准；
2. 分数以三个任务的当前 `results_new/score_rerun_results.txt` 为准；
3. 带 `before_deepvelo_sim_rerun_20260613` 的目录明确属于旧版；
4. TIvelo `filter_genes` 目录属于后续诊断，不是 Figure 3 输入。

### Real branching/MURK genes

1. 当前结果位于 `plot/final_genes/`；
2. `before_gene_cbdir_rerun_20260613_013621/final_genes/` 是旧版；
3. 对最终 PDF 有疑问时，应比较图中显示分数与两个 `score.txt`，而不是比较目录时间戳。

### Kinetics analysis

1. 当前 release 统一结果位于 `evaluation/kinetics_recovery/results/final_seed123/`；
2. DeepKINET 的正式结果使用 seed 123，参数恢复只读取 `degradation_rate` 作为 gamma；
3. Figure 5B–F 的可复现候选素材位于 `evaluation/kinetics_recovery/figures/final_seed123/`；
4. Figure 5A 仍是人工选择的 RPE1 流线素材，暂不纳入本次 seed123 更新。

### Ensemble

1. 最终方法为 row-wise geometric median，论文简称 GeoMedian-4；
2. Mean-4 只作为消融基线；
3. S10 直接读取 `scores.csv`，S12 直接读取四张 stream PNG；S11 的当前论文构建脚本仍读取带 `with_row_geometric_median` 历史内部命名的两张 PNG；
4. 发布仓库中统一图件命名时，必须同时修改构建脚本并进行 PDF 视觉核验；
5. 重跑或引用前先确认 Git 工作区和具体 commit。

## 8. 方法与数据集命名契约

论文正文、图、图注和表格中的方法名与数据集名，统一以 `manifests/naming_conventions.md` 为唯一规范。运行脚本、结果路径和历史记录继续使用 operational ID 或旧别名，不据此修改文件路径。

特别注意：`rpe1`（Battich scEU-seq）与 `rpe1_gen`（展示名 `rpe1 (VeloCycle)`）是两个不同数据集，不能合并或互作别名。绘图脚本如需本地映射，必须与该规范逐项一致，不得引入冲突别名。

## 9. 已知 provenance 风险

1. `/mnt/second19T/zhaozelin/velocity` 整体不是 Git 仓库，只有 `transition_ensemble/` 具有独立 Git 历史；
2. 少数当前脚本的修改时间晚于最终结果，无法仅凭 mtime 恢复精确代码快照；
3. Figure 1 和 Figure 6 缺少完整的可编辑或自动排版源；
4. Figure 5A 的有效素材位于名称上看似可删除的 archive 目录；
5. Figure 5、Figure 6 和部分补充图的最终 PDF 包含人工组合、排版或裁剪步骤；S9–S12 已有统一论文侧构建脚本；
6. 某些目录同时保留 canonical、legacy reproduction、diagnostic 和 before-rerun 文件，名称相似但作用不同。

因此，不应仅凭“文件较旧”“目录名含 old/archive/before”或“已有 Git”决定删除文件。

## 10. 新 agent 的推荐阅读顺序

开始新任务时建议依次执行：

1. 阅读本文件，确定任务属于 cell-level、simulation、real gene-level、kinetic 或 ensemble；
2. 以 Overleaf `merged.tex` 和最终 PDF 为论文事实来源；
3. 从最终 panel 反向定位到直接素材，再定位绘图脚本和指标输入；
4. 查找相邻的 `manifest.json`、`render_manifest.csv`、`scores.csv` 或 `score.txt`；
5. 对存在当前/旧版并列的结果，比较图中可见数值或文件哈希；
6. 如果涉及 `transition_ensemble`，先检查 Git 状态；
7. 如果任务只是排查或解释，不要重跑实验、修改文件或清理历史目录；
8. 修改前记录输入、输出、环境和版本；修改后再验证最终 PDF，而不是只看脚本运行成功。

推荐的只读入口命令：

```bash
ssh 122.207.108.8
cd /mnt/second19T/zhaozelin/velocity

find indicators plot transition_ensemble -maxdepth 2 -type d | sort
find plot -type d \( -iname '*old*' -o -iname '*archive*' -o -iname '*before*' \) | sort

cd transition_ensemble
git status --short
git log -5 --oneline
```

## 11. 清理文件前的硬性约束

在建立完整的 final-asset manifest 之前：

- 不要整体删除 `plot/archive_old_versions_20260613/`；
- 不要删除 Figure 5A 的五张 RPE1 stream 图；
- 不要用当前脚本覆盖已有结果来“验证”版本；
- 不要把 `results_new_before_*` 重新命名成当前结果；
- 不要清理 ensemble 的旧命名图件，除非已同步修改 S11 构建脚本并核对 Overleaf 最终版本；
- 不要假设人工拼版图可以由单个 Python 脚本完整重建。

更安全的后续整理方式是先建立一份包含以下字段的 manifest：

```text
figure
panel
final_asset
generator_script
direct_inputs
upstream_results
generated_at
sha256
git_commit_or_unversioned
status_current_or_historical
notes
```

只有当所有最终 panel 都有明确的 canonical asset 后，才适合继续迁移或删除历史目录。
