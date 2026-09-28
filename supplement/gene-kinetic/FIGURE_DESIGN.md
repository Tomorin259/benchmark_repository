# Gene-kinetics phase-portrait design

本目录最终批量只交付预定基因和边的 masked gene spliced–unspliced phase portraits。每个组合输出九方法 2×5 面板：前 9 格为方法，第 10 格为共享 cluster 颜色图例。方法名下显示预定边的完整名称和逐边 gene-CBDir 分值；缺失分数显示 `完整边名: N/A`。

相图调用 `evaluation/gene/plotting/phase_portrait.py:266`，按目标 cluster 和单基因切片后传入。固定参数为 `Ms/Mu`、`velocity/velocity_u`、`color_key="clusters"`、`arrow_grid=(25,25)`、`s=12`、`alpha=0.6`、默认 `reference_length=0.05` 和 `max_length=0.10`。坐标轴和刻度隐藏。显示参数不改变速度层、输入数据或 CBDir。

UMAP 与 velocity stream 不属于 gene-level 成品；它们使用整体细胞速度场，目标 cluster 仅用于构造 masked AnnData。

## 完整批量清单

| 数据集 | 基因 | mask cluster | 预定 score_edges |
| --- | --- | --- | --- |
| Gastrulation | Hba-x, Gypa, Alas2, Fam210b, Klf3, Hemgn, Fech | Erythroid1, Erythroid2, Erythroid3 | Erythroid1→Erythroid2; Erythroid2→Erythroid3 |
| DentateGyrus_lamanno | Ncald | Nbl2, ImmGranule1, ImmGranule2, Granule | ImmGranule1→ImmGranule2; ImmGranule2→Granule |
| DentateGyrus_lamanno | Diaph3 | nIPC, Nbl1, Nbl2 | nIPC→Nbl1; Nbl1→Nbl2 |
| DentateGyrus_lamanno | Dab1 | CA, CA1-Sub, CA2-3-4 | CA→CA1-Sub |
| capo_pre | E2F7, FAM111A | G1, S, G2M | G1→S; S→G2M |
| retina | Celf4 | Neuroblast, RGC, AC/HC | Neuroblast→RGC; Neuroblast→AC/HC |
| mousebrain | Flrt2 | V-SVZ, Upper Layer, Deeper Layer | V-SVZ→Upper Layer |
| pancreas | G6pc2 | Pre-endocrine, Beta, Alpha | Pre-endocrine→Beta; Pre-endocrine→Alpha |

Diaph3 的 `nIPC→Nbl1` 使用本目录补算 CSV；DeepVelo_scv 的 gene-missing 保持 N/A。其他分值直接读取 evaluation 逐边 CSV。每个组合保存 PNG、PDF、SVG 和 manifest，记录输入 locator、mask、方法、边分数、预处理诊断和显示参数。
