# Gene kinetics and local cell-direction study

本目录存放针对少数有生物学先验的基因与发育边的补充分析。研究问题是：方法在特定基因、特定阶段恢复的 RNA 动力学方向，是否与同一边上的细胞方向表现相联系；哪些例外提示模型假设、基因尺度或投影流程仍有问题。本目录不替代全数据集排名，也不更改论文或现有主评估结果。

## 数据与判定边界

- 九个 gene-analysis 方法的结果与有效性以 `../../evaluation/gene/results/gene_cbdir_reference30_all.csv` 为准；逐边值见同目录的 `gene_cbdir_edge_details_reference30_all.csv`。只将 `status=ok` 且分数有限的方法计入覆盖数。缺失方法显示为缺失，不记为零分。
- 固定 30 邻居、cluster 标签、velocity 层及边定义沿用 `../../evaluation/gene/run_reference30_gene_cbdir.py`。方法原生的 `Ms`、`Mu`、spliced/unspliced velocity 保持不变。图中表达趋势只是当前数据的描述，不作为独立生物学先验。
- 细胞层面的逐边 `highdim_cbdir` 与 `transcosine` 已在 `../../evaluation/cell/results/cell_metrics_edge_details.csv`；DeepVelo、UniTVelo、veloVI 的 gene-output-matched 逐边结果在 `../../evaluation/cell/gene_matched/results/cell_metrics_edge_details.csv`。比较之前必须逐方法核对源 H5AD、velocity key、方法配置与边的源细胞集合；不能仅按方法显示名拼接。
- 选边先看原始研究及功能实验，再看当前数据的表达、邻接和可评估细胞数。逐边分数只能检验预定假设，不能反向决定哪些边是“生物学相关边”。

## 预定的基因—边分析

| 优先级 | 数据集与基因 | 研究边 | 对照或限制 | 先验 |
| --- | --- | --- | --- | --- |
| 核心：晚期动力学切换 | Gastrulation：Hba-x (7/9)、Gypa (8/9)、Alas2 (7/9)、Fam210b (8/9)，并列 Klf3、Hemgn、Fech (各 9/9) | Erythroid2 → Erythroid3 | Erythroid1 → Erythroid2 为前期对照；不将 Blood progenitor 边解释为同一 MURK 事件 | Barile 等以 Erythroid3 相对 Erythroid2 的 unspliced–spliced 斜率升高定义 mouse MURK；Hba-x 是原文相图例子，Fam210b 在 Gata1 嵌合实验中展示晚期变化。以上基因均在原研究 MURK 清单中。[原文](https://doi.org/10.1186/s13059-021-02414-y) |
| 核心：颗粒细胞成熟 | DentateGyrus_lamanno：Ncald (9/9) | ImmGranule1 → ImmGranule2；ImmGranule2 → Granule | CA → CA1-Sub 与 CA → CA2-3-4 仅作区域/分支补充比较 | 原 cellDancer 研究将 Ncald 列为分支动力学基因；敲除研究关联齿状回神经发生，并指出 DG、CA3 表达较丰富。[cellDancer](https://doi.org/10.1038/s41587-023-01728-5)；[功能研究](https://doi.org/10.3389/fnmol.2019.00019) |
| 核心：祖细胞退出 | DentateGyrus_lamanno：Diaph3 (8/9) | nIPC → Nbl1；Nbl1 → Nbl2 为次级检查 | **nIPC → Nbl1 尚未进入现有 gene 逐边表，必须单独补算。**固定邻居图中 304 个 nIPC 源细胞邻接 Nbl1。 | 原 cellDancer 研究列为分支基因；敲除实验主要证明神经祖细胞分裂/神经发生作用，研究组织是皮层，不能直接声称已在齿状回做因果验证。[cellDancer](https://doi.org/10.1038/s41587-023-01728-5)；[敲除研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC8102060/)；[TIvelo 的该数据集边定义](https://tivelo.readthedocs.io/en/latest/notebooks/notebooks/DentateGyrus2.html) |
| 核心：细胞周期机制对照 | capo_pre：E2F7、FAM111A (均 9/9) | G1 → S；S → G2M | G2M → G1 可作闭环检查，不作为复制起始的主要边 | E2F7 调控周期基因与 S 期；FAM111A 促进复制起始。[E2F7](https://pmc.ncbi.nlm.nih.gov/articles/PMC3333892/)；[FAM111A](https://pmc.ncbi.nlm.nih.gov/articles/PMC10551639/) |
| 探索：迁移/定位 | DentateGyrus_lamanno：Dab1 (8/9) | CA → CA1-Sub | CA1 迁移表型不等于 Dab1 转录本的特殊动力学 | [Dab1 海马实验](https://pmc.ncbi.nlm.nih.gov/articles/PMC10089055/) |
| 探索：视网膜分化表达 | retina：Celf4 (9/9) | Neuroblast → RGC；Neuroblast → AC/HC 同时展示 | 原位表达证据不能证明命运决定或复杂动力学；两边各仅有 52、41 个可评估源细胞 | [Celf4 时空表达](https://pmc.ncbi.nlm.nih.gov/articles/PMC3825781/) |
| 探索：皮层迁移 | mousebrain：Flrt2 (9/9) | V-SVZ → Upper Layer | FLRT2 是皮层板提供的迁移信号；它自己的 RNA velocity 不是迁移作用的直接读数。Upper Layer → Deeper Layer 暂不作发育方向解释。 | [FLRT2/Unc5D 实验](https://doi.org/10.1038/emboj.2011.189) |
| 探索：终末功能 | pancreas：G6pc2 (9/9) | Pre-endocrine → Beta；→ Alpha 为对照 | β 细胞功能证据不等于 β 细胞命运决定证据 | [G6pc2 敲除](https://pmc.ncbi.nlm.nih.gov/articles/PMC7331801/) |

Snap25、Hcn1、Igf1r、LMNA 暂不进入主要边分析：现有文献主要支持较晚功能或不同过程，未充分锚定本评估图中的过渡。若作为探索性图，分别只考虑 retina 的 Neuroblast → PR、pancreas 的 Pre-endocrine → Beta、bonemarrow 的 Ery_1 → Ery_2，并在图注明确证据局限。其他 cluster 边不因已有数值可用就自动纳入。

## 执行顺序

1. **冻结分析清单。**逐基因记录有效方法、缺失方法、对应 H5AD 与速度层、所选边和可评估源细胞数。主分析与探索分析分开；边的选取不随分数变动。
2. **检查动力学形态。**在每条预定边上展示同一批细胞的 spliced/unspliced 相图、方法速度箭头和分群表达。红系同时展示 Erythroid1/2/3，以核查晚期斜率变化；齿状回沿分支展示，避免全图混合造成假象。图仅用仓库已有的可审计相图绘制约定。
3. **比较局部方向。**从已有逐边表提取 gene-CBDir 与 cell-level 的 `highdim_cbdir`、`transcosine`；对照覆盖数和源细胞数。它们衡量不同空间中的方向，须分别报告，不合成一个“正确性”分数。Diaph3 的 nIPC → Nbl1 使用相同固定邻居和原指标公式补算到本目录，不回填或覆盖原始表。
4. **检验跨层联系。**在相同方法输出、相同数据集和边上比较基因方向与细胞方向；同时展示不一致的方法。相关性只能支持关联，不能证明一个基因导致整群细胞方向好坏。若红系结果提示 MURK 敏感性，再预定使用同一邻居与投影设置做 MURK 基因集移除实验，并把它标为敏感性分析，而非模型修复。
5. **独立出图与核验。**每个生物学案例输出独立子图与源数据表；图注明示文献先验、边、细胞数、有效方法数及缺失方法。核对相图方向、数据来源和结果表一致后，再决定是否组图或写入论文。所有论文文字先形成建议稿，不在本研究步骤中直接改在线稿。

## 结果存放约定

本轮试画最终只交付基因层面的 masked spliced--unspliced 相图。UMAP/velocity stream
使用的是整体细胞速度场，不能作为单基因动力学结果，因此不纳入本目录的成品图；目标
cluster 仅用于构造 masked AnnData 后送入仓库的 `evaluation/gene/plotting/phase_portrait.py`。

后续生成的紧凑逐边表放在本目录 `results/`，独立图放在 `figures/`，每个结果附输入定位、参数、方法覆盖和 SHA256 或等效来源记录。只在实际产出时创建这些子目录；不复制现有评估模块的整套表或大型 H5AD。论文主结果与补充研究结果保持分开。

首轮相图恢复仓库默认显示参数 `reference_length=0.05`、`max_length=0.10`，并统一使用 `s=12`、`alpha=0.6`；仅影响显示，不改变基因速度层、细胞数据或任何 CBDir 分数。
