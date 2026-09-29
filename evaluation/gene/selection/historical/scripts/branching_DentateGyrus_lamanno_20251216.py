import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
from scipy import stats
import warnings

warnings.filterwarnings('ignore')

# ==============================================================================
# 1. MURK 特有的“双重基因筛选”预处理
# ==============================================================================
print("[Step 1] 执行双重基因筛选 (Scanpy HVG + scVelo Genes)...")
adata=sc.read('velocity/data/DentateGyrus_lamanno/DentateGyrus_lamanno_for_scvelo.h5ad')
# 1.0 备份原始数据 (因为 Scanpy 和 scVelo 的处理会冲突，必须分开跑)
# 假设 adata 是刚读取进来的 raw counts
adata_raw = adata.copy()

# --- A. Scanpy 路线：找高变基因 (HVG) ---
print("   -> Running Scanpy HVG selection...")
# 这一步是为了找“表达量差异大”的基因
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.05, min_disp=0.1)
list_hvg = adata.var_names[adata.var.highly_variable].tolist()

# --- B. scVelo 路线：找动力学基因 ---
print("   -> Running scVelo gene selection...")
# 必须回滚到原始数据！否则会重复归一化
adata = adata_raw.copy() 
# scVelo 的这个函数包含了 normalize 和 log1p，同时筛选适合 Velocity 的基因
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
list_scv = adata.var_names.tolist()

# --- C. 取并集 (The MURK Union) ---
# 既要高变，也要有动力学信号
# set() 去重，然后转回 list
final_gene_list = list(set(list_hvg + list_scv))
print(f"   -> Scanpy HVG: {len(list_hvg)}, scVelo Genes: {len(list_scv)}")
print(f"   -> 合并后总基因数: {len(final_gene_list)}")

# --- D. 最终数据准备 ---
# 再次回滚，只保留并集基因，然后做最终的平滑
adata = adata_raw.copy()
# 这一步是必须的：因为刚才只是为了拿名单，现在要对选定的基因做真正的预处理
# 注意：这里不再筛选 n_top_genes，而是指定我们要的基因列表
adata = adata[:, final_gene_list].copy()

# 执行 scVelo 标准化 (归一化 + 对数化)
scv.pp.filter_and_normalize(adata, min_shared_counts=20) 
# 计算 Moments (平滑) -> 生成 Ms/Mu
scv.pp.moments(adata, n_neighbors=30, n_pcs=30)

# ==============================================================================
# 2. 准备计算矩阵
# ==============================================================================
# 此时 adata.var_names 已经是我们精心挑选的 final_gene_list
valid_genes = adata.var_names

if hasattr(adata.layers['Ms'], 'toarray'):
    df_Ms = pd.DataFrame(adata.layers['Ms'].toarray(), index=adata.obs_names, columns=valid_genes)
    df_Mu = pd.DataFrame(adata.layers['Mu'].toarray(), index=adata.obs_names, columns=valid_genes)
else:
    df_Ms = pd.DataFrame(adata.layers['Ms'], index=adata.obs_names, columns=valid_genes)
    df_Mu = pd.DataFrame(adata.layers['Mu'], index=adata.obs_names, columns=valid_genes)

# ==============================================================================
# 3. MURK 核心统计函数 (保持不变)
# ==============================================================================
def calculate_murk_statistics(mask, global_max_s, global_max_u, confidence=0.95):
    n = np.sum(mask)
    if n < 10: return None

    # 提取与 Scaling
    X = df_Ms.loc[mask].values
    Y = df_Mu.loc[mask].values
    
    s_denom = global_max_s.values.copy(); s_denom[s_denom == 0] = 1.0
    u_denom = global_max_u.values.copy(); u_denom[u_denom == 0] = 1.0
    
    x = X / s_denom
    y = Y / u_denom

    # OLS 回归
    mx = np.mean(x, axis=0)
    my = np.mean(y, axis=0)
    x_c = x - mx
    y_c = y - my

    numerator = np.sum(x_c * y_c, axis=0)
    denominator = np.sum(x_c ** 2, axis=0)

    with np.errstate(divide='ignore', invalid='ignore'):
        slope = numerator / denominator
        slope[denominator == 0] = 0.0

    # 误差与区间
    intercept = my - slope * mx
    y_pred = slope * x + intercept
    sse = np.sum((y_pred - y) ** 2, axis=0)

    with np.errstate(divide='ignore', invalid='ignore'):
        se = np.sqrt(sse / (n - 2) / denominator)
        se[denominator == 0] = np.inf

    t_val = stats.t.ppf((1 + confidence) / 2, n - 2)
    
    return {
        'slope': slope,
        'k_min': slope - t_val * se,
        'k_max': slope + t_val * se
    }

# ==============================================================================
# 4. 执行分析任务 (保持之前的逻辑)
# ==============================================================================

# --- Task 1: Neuron vs Glia ---
print("\nRunning Task 1: Neuron vs Glia...")
cloud_neuron = ['RadialGlia', 'nIPC', 'Nbl1', 'Nbl2', 'ImmGranule1', 'ImmGranule2', 'Granule', 'CA', 'CA1-Sub', 'CA2-3-4']
cloud_glia   = ['RadialGlia', 'GlialProg', 'ImmAstro', 'OPC']

mask_n = adata.obs['clusters'].isin(cloud_neuron)
mask_g = adata.obs['clusters'].isin(cloud_glia)
mask_union = mask_n | mask_g

max_s = df_Ms.loc[mask_union].max(axis=0)
max_u = df_Mu.loc[mask_union].max(axis=0)

stats_n = calculate_murk_statistics(mask_n, max_s, max_u)
stats_g = calculate_murk_statistics(mask_g, max_s, max_u)

if stats_n and stats_g:
    gap_n = np.maximum(stats_n['k_min'] - stats_g['k_max'], 0)
    gap_g = np.maximum(stats_g['k_min'] - stats_n['k_max'], 0)
    
    df_res1 = pd.DataFrame({
        'gene': valid_genes,
        'separation_gap': gap_n + gap_g,
        'slope_neuron': stats_n['slope'],
        'slope_glia': stats_g['slope']
    }).sort_values('separation_gap', ascending=False)
    
    df_res1.to_csv('velocity/indicators/gene/DentateGyrus_lamanno_analyze/branching_genes/Result_Neuron_vs_Glia.csv', index=False)
    print(f"Task 1 Top: {df_res1.iloc[0]['gene']}")

# --- Task 2: Granule vs Pyramidal ---
print("\nRunning Task 2: Granule vs Pyramidal...")
cloud_granule   = ['Nbl2', 'ImmGranule1', 'ImmGranule2', 'Granule']
cloud_pyramidal = ['Nbl2', 'CA', 'CA1-Sub', 'CA2-3-4']

mask_gr = adata.obs['clusters'].isin(cloud_granule)
mask_py = adata.obs['clusters'].isin(cloud_pyramidal)
mask_union_2 = mask_gr | mask_py

max_s_2 = df_Ms.loc[mask_union_2].max(axis=0)
max_u_2 = df_Mu.loc[mask_union_2].max(axis=0)

stats_gr = calculate_murk_statistics(mask_gr, max_s_2, max_u_2)
stats_py = calculate_murk_statistics(mask_py, max_s_2, max_u_2)

if stats_gr and stats_py:
    gap_gr = np.maximum(stats_gr['k_min'] - stats_py['k_max'], 0)
    gap_py = np.maximum(stats_py['k_min'] - stats_gr['k_max'], 0)
    
    df_res2 = pd.DataFrame({
        'gene': valid_genes,
        'separation_gap': gap_gr + gap_py,
        'slope_granule': stats_gr['slope'],
        'slope_pyramidal': stats_py['slope']
    }).sort_values('separation_gap', ascending=False)
    
    df_res2.to_csv('velocity/indicators/gene/DentateGyrus_lamanno_analyze/branching_genes/Result_Granule_vs_Pyramidal.csv', index=False)
    print(f"Task 2 Top: {df_res2.iloc[0]['gene']}")