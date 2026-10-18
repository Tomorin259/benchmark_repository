# Gene kinetics and local cell-direction cases

本目录保存少量预先定义的基因—发育边案例及跨层描述性比较。
完整数据集的结果见 [Cell](../../evaluation/cell/README.md) 和
[Gene](../../evaluation/gene/README.md)。

## 已有材料与运行入口

| 位置 | 内容 |
| --- | --- |
| [figures/](figures/) | masked spliced–unspliced 相图、逐案例 manifest |
| [analysis/](analysis/) | 配对来源、逐边观测、聚合、相关性和方向一致性表 |
| [analysis/source_mapping_audit.csv](analysis/source_mapping_audit.csv) | 方法 H5AD、速度层和案例来源对应 |
| [results/](results/) | 独立 Cell 配对结果及 Diaph3 额外边的补算 |
| [相图显示与切片](#相图显示与切片) | 展示边、cluster 掩膜和显示参数 |

已有紧凑材料可用于重新汇总；从仓库根目录输出到新目录：

```bash
python supplement/gene-kinetic/analysis/gene_kinetic_analysis.py \
  --repo-root . --output-dir outputs/gene-edge-analysis
```

重新绘制相图另需方法 H5AD、参考标签和科学绘图环境。
[gene_kinetic_pilot.py](gene_kinetic_pilot.py) 要求显式指定 `--output-dir`，
将相图和每例 manifest 写入仓库外或 `outputs/` 下的空目录；原图件与输入不改写。

```bash
python supplement/gene-kinetic/gene_kinetic_pilot.py \
  --artifact-root /path/to/velocity-artifacts \
  --case dentate_ncald --output-dir outputs/gene-cases
```

整群 UMAP 和 velocity stream 属于细胞速度场，不作为单基因相图。

## 比较口径

Gene CBDir 来自 [reference-30 逐边表](../../evaluation/gene/results/gene_cbdir_edge_details_reference30_all.csv)；
Cell 比较使用已对齐的对应方法输出，DeepVelo、UniTVelo、veloVI 的
Gene-output-matched 来源与主 Cell 来源分开。有效方法要求 `status=ok`
且分数有限；缺失值不记零。

比较前核对 H5AD、速度层、方法配置和源细胞集合。同一边上的多个基因
先按 edge–method 聚合，避免将共享细胞方向的基因当作独立重复。
`batc_candidate_source_mean` 是分支最大化前的候选均值，
`branch_winner_share` 是候选获胜比例，两者均不同于最终 BATC。

这些案例比较预定发育边上的 Gene–Cell 分数关联和方向一致性。
[分析报告](analysis/analysis_report_20260920.md) 汇总逐边结果、方法覆盖和聚合单位。

## 案例与文献依据

下表保留案例定义、先验来源和解释范围。选择依据不由当前分数反推。

| 分析范围 | 数据集与基因 | 研究边 | 对照或限制 | 先验 |
| --- | --- | --- | --- | --- |
| 核心：晚期动力学切换 | Gastrulation：Hba-x (7/9)、Gypa (8/9)、Alas2 (7/9)、Fam210b (8/9)，并列 Klf3、Hemgn、Fech (各 9/9) | Erythroid2 → Erythroid3 | Erythroid1 → Erythroid2 为前期对照；不将 Blood progenitor 边解释为同一 MURK 事件 | Barile 等以 Erythroid3 相对 Erythroid2 的 unspliced–spliced 斜率升高定义 mouse MURK；Hba-x 是原文相图例子，Fam210b 在 Gata1 嵌合实验中展示晚期变化。以上基因均在原研究 MURK 清单中。[原文](https://doi.org/10.1186/s13059-021-02414-y) |
| 核心：颗粒细胞成熟 | DentateGyrus_lamanno：Ncald (9/9) | ImmGranule1 → ImmGranule2；ImmGranule2 → Granule | CA → CA1-Sub 与 CA → CA2-3-4 仅作区域/分支补充比较 | 原 cellDancer 研究将 Ncald 列为分支动力学基因；敲除研究关联齿状回神经发生，并指出 DG、CA3 表达较丰富。[cellDancer](https://doi.org/10.1038/s41587-023-01728-5)；[功能研究](https://doi.org/10.3389/fnmol.2019.00019) |
| 核心：祖细胞退出 | DentateGyrus_lamanno：Diaph3 (8/9) | nIPC → Nbl1；Nbl1 → Nbl2 为次级检查 | nIPC → Nbl1 的独立补算见本目录 `results/diaph3_nipc_to_nbl1.csv`，主评估表保持原样。固定邻居图中 304 个 nIPC 源细胞邻接 Nbl1。 | 原 cellDancer 研究列为分支基因；敲除实验主要证明神经祖细胞分裂/神经发生作用，研究组织为皮层。[cellDancer](https://doi.org/10.1038/s41587-023-01728-5)；[敲除研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC8102060/)；[TIvelo 的该数据集边定义](https://tivelo.readthedocs.io/en/latest/notebooks/notebooks/DentateGyrus2.html) |
| 核心：细胞周期机制对照 | capo_pre：E2F7、FAM111A (均 9/9) | G1 → S；S → G2M | G2M → G1 可作闭环检查，不作为复制起始的主要边 | E2F7 调控周期基因与 S 期；FAM111A 促进复制起始。[E2F7](https://pmc.ncbi.nlm.nih.gov/articles/PMC3333892/)；[FAM111A](https://pmc.ncbi.nlm.nih.gov/articles/PMC10551639/) |
| 探索：迁移/定位 | DentateGyrus_lamanno：Dab1 (8/9) | CA → CA1-Sub | CA1 迁移表型不等于 Dab1 转录本的特殊动力学 | [Dab1 海马实验](https://pmc.ncbi.nlm.nih.gov/articles/PMC10089055/) |
| 探索：视网膜分化表达 | retina：Celf4 (9/9) | Neuroblast → RGC；Neuroblast → AC/HC 同时展示 | 原位研究提供分化阶段的表达证据；两边分别有 52、41 个可评估源细胞 | [Celf4 时空表达](https://pmc.ncbi.nlm.nih.gov/articles/PMC3825781/) |
| 探索：皮层迁移 | mousebrain：Flrt2 (9/9) | V-SVZ → Upper Layer | FLRT2 是皮层板提供的迁移信号；它自己的 RNA velocity 不是迁移作用的直接读数。Upper Layer → Deeper Layer 暂不作发育方向解释。 | [FLRT2/Unc5D 实验](https://doi.org/10.1038/emboj.2011.189) |
| 探索：终末功能 | pancreas：G6pc2 (9/9) | Pre-endocrine → Beta；→ Alpha 为对照 | 敲除研究提供 β 细胞功能证据 | [G6pc2 敲除](https://pmc.ncbi.nlm.nih.gov/articles/PMC7331801/) |

Snap25、Hcn1、Igf1r、LMNA 的文献证据未被充分锚定到此案例集合的预定过渡，
因此不列为这里的主要边分析；不影响它们在主 Gene 评估中的既有用途。
## 相图显示与切片

上表定义分析问题和先验；下表列出同组全部基因的显示掩膜与标注分数的边。
显示中同时纳入表列的对照或分支边。

| 案例组（对应上表） | mask clusters | 展示分数的边 |
| --- | --- | --- |
| 晚期动力学切换 | Erythroid1, Erythroid2, Erythroid3 | Erythroid1→Erythroid2; Erythroid2→Erythroid3 |
| 颗粒细胞成熟 | Nbl2, ImmGranule1, ImmGranule2, Granule | ImmGranule1→ImmGranule2; ImmGranule2→Granule |
| 祖细胞退出 | nIPC, Nbl1, Nbl2 | nIPC→Nbl1; Nbl1→Nbl2 |
| 细胞周期机制对照 | G1, S, G2M | G1→S; S→G2M |
| 迁移/定位 | CA, CA1-Sub, CA2-3-4 | CA→CA1-Sub |
| 视网膜分化表达 | Neuroblast, RGC, AC/HC | Neuroblast→RGC; Neuroblast→AC/HC |
| 皮层迁移 | V-SVZ, Upper Layer, Deeper Layer | V-SVZ→Upper Layer |
| 终末功能 | Pre-endocrine, Beta, Alpha | Pre-endocrine→Beta; Pre-endocrine→Alpha |

每个案例为九方法 2×5 面板，第十格为共享 cluster 图例；方法名下标注完整边名
和逐边 Gene CBDir，缺失分数为 `N/A`。Diaph3 的 nIPC→Nbl1 使用独立补算
[CSV](results/diaph3_nipc_to_nbl1.csv)；DeepVelo_scv 缺失该基因时保持 N/A，
其他边读取主评估逐边表。

相图经 [plot_phase_portrait](../../evaluation/gene/plotting/phase_portrait.py)
绘制，保留原生 `Ms/Mu`，将方法对应的速度分量映射为
`velocity/velocity_u`。显示参数为 `color_key="clusters"`、
`arrow_grid=(25,25)`、log1p 箭长、参考分位数 0.90、
`reference_length=0.05`、`max_length=0.10`、`s=12`、`alpha=0.6`；
坐标轴和刻度隐藏。这些显示设置不改变源 H5AD、速度层或 CBDir。
每例保存 PNG、PDF 和 manifest，记录输入 locator、掩膜、方法、边分数、
预处理诊断及显示参数。
