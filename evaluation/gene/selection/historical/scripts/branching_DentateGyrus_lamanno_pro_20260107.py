import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
from scipy import stats
import warnings

warnings.filterwarnings('ignore')

# ==============================================================================
# 1. 预处理：双重基因筛选 (保持原方案)
# ==============================================================================
print("[Step 1] 执行双重基因筛选 (Scanpy HVG + scVelo Genes)...")
# 读取数据
adata = sc.read('velocity/data/DentateGyrus_lamanno/DentateGyrus_lamanno.h5ad')
adata_raw = adata.copy()

# A. Scanpy HVG
print("   -> Running Scanpy HVG selection...")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.05, min_disp=0.1)
list_hvg = adata.var_names[adata.var.highly_variable].tolist()

# B. scVelo Genes
print("   -> Running scVelo gene selection...")
adata = adata_raw.copy()
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
list_scv = adata.var_names.tolist()

# C. 取并集
final_gene_list = list(set(list_hvg + list_scv))
print(f"   -> 合并后总基因数: {len(final_gene_list)}")

# D. 最终数据准备
adata = adata_raw.copy()
adata = adata[:, final_gene_list].copy()
scv.pp.filter_and_normalize(adata, min_shared_counts=20)
scv.pp.moments(adata, n_neighbors=30, n_pcs=30)

# ==============================================================================
# 2. 准备计算矩阵
# ==============================================================================
valid_genes = adata.var_names
if hasattr(adata.layers['Ms'], 'toarray'):
    df_Ms = pd.DataFrame(adata.layers['Ms'].toarray(), index=adata.obs_names, columns=valid_genes)
    df_Mu = pd.DataFrame(adata.layers['Mu'].toarray(), index=adata.obs_names, columns=valid_genes)
else:
    df_Ms = pd.DataFrame(adata.layers['Ms'], index=adata.obs_names, columns=valid_genes)
    df_Mu = pd.DataFrame(adata.layers['Mu'], index=adata.obs_names, columns=valid_genes)

# ==============================================================================
# 3. 核心函数定义
# ==============================================================================

def calculate_branching_slope(mask, global_max_s, global_max_u, confidence=0.95):
    """
    计算相图斜率 (k) 及其置信区间 (k_min, k_max)
    """
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

def calculate_trend_spearman(mask, rank_dict, obs_clusters):
    """
    计算表达量趋势: Spearman Correlation (Ms vs Rank)
    """
    n = np.sum(mask)
    if n < 10: return np.zeros(len(valid_genes))

    current_cells = obs_clusters[mask]
    
    # 映射 Rank
    # 填充可能出现的 NaN (如果有未定义的簇) 为 0
    x_ranks = current_cells.map(rank_dict).fillna(0).values
    Y_expr = df_Ms.loc[mask].values

    rhos = []
    for i in range(Y_expr.shape[1]):
        y = Y_expr[:, i]
        if np.std(y) == 0:
            rhos.append(0.0)
        else:
            rho, _ = stats.spearmanr(x_ranks, y)
            if np.isnan(rho): rho = 0.0
            rhos.append(rho)
            
    return np.array(rhos)

# ==============================================================================
# 4. 执行分析任务：Granule vs Pyramidal
# ==============================================================================
print("\nRunning Task 2: Granule vs Pyramidal...")

# 定义参数
SPEARMAN_THRESHOLD = 0.1  # [关键参数] 趋势相关性过滤阈值

cloud_granule   = ['Nbl2', 'ImmGranule1', 'ImmGranule2', 'Granule']
cloud_pyramidal = ['Nbl2', 'CA', 'CA1-Sub', 'CA2-3-4']

# Rank 定义
rank_dict_granule = {
    'Nbl2': 1, 
    'ImmGranule1': 2, 
    'ImmGranule2': 3, 
    'Granule': 4
}
rank_dict_pyramidal = {
    'Nbl2': 1, 
    'CA': 2, 
    'CA1-Sub': 3, 
    'CA2-3-4': 3
}

# 1. Mask 准备
mask_gr = adata.obs['clusters'].isin(cloud_granule)
mask_py = adata.obs['clusters'].isin(cloud_pyramidal)
mask_union = mask_gr | mask_py

# 2. 全局缩放因子
max_s = df_Ms.loc[mask_union].max(axis=0)
max_u = df_Mu.loc[mask_union].max(axis=0)

# 3. 计算动力学 (Slope)
print("   -> Calculating Slopes...")
stats_gr = calculate_branching_slope(mask_gr, max_s, max_u)
stats_py = calculate_branching_slope(mask_py, max_s, max_u)

# 4. 计算趋势 (Spearman)
print("   -> Calculating Spearman Trends...")
rho_gr = calculate_trend_spearman(mask_gr, rank_dict_granule, adata.obs['clusters'])
rho_py = calculate_trend_spearman(mask_py, rank_dict_pyramidal, adata.obs['clusters'])

# 5. 整合与筛选 (Gatekeeper Strategy)
if stats_gr and stats_py:
    # A. 计算动力学差异 (Slope Gap)
    gap_gr = np.maximum(stats_gr['k_min'] - stats_py['k_max'], 0)
    gap_py = np.maximum(stats_py['k_min'] - stats_gr['k_max'], 0)
    total_slope_gap = gap_gr + gap_py
    
    # B. 计算趋势置信度 (Max Spearman)
    trend_validity = np.maximum(np.abs(rho_gr), np.abs(rho_py))
    
    # C. 创建临时 DataFrame
    df_temp = pd.DataFrame({
        'gene': valid_genes,
        'slope_gap': total_slope_gap,
        'trend_validity': trend_validity,
        'slope_granule': stats_gr['slope'],
        'slope_pyramidal': stats_py['slope'],
        'rho_granule': rho_gr,
        'rho_pyramidal': rho_py
    })
    
    # D. 执行过滤 (Gatekeeper)
    # 只有至少在一个分支中有显著趋势的基因才被保留
    print(f"   -> Filtering genes with Trend Validity < {SPEARMAN_THRESHOLD}...")
    df_filtered = df_temp[df_temp['trend_validity'] > SPEARMAN_THRESHOLD].copy()
    
    # E. 执行排序
    # 在保留下来的基因中，纯粹按动力学差异大小排序
    df_res2 = df_filtered.sort_values('slope_gap', ascending=False)
    
    # 6. 保存
    out_path = 'velocity/indicators/gene/DentateGyrus_lamanno_analyze/branching_genes/select_results_filtered.csv'
    df_res2.to_csv(out_path, index=False)
    
    print(f"   -> 原始基因数: {len(df_temp)}, 过滤后基因数: {len(df_res2)}")
    if not df_res2.empty:
        print(f"Task 2 Top Gene: {df_res2.iloc[0]['gene']}")
    print(f"Result saved to: {out_path}")