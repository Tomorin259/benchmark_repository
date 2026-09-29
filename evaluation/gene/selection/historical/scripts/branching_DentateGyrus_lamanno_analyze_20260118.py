import scanpy as sc
import scvelo as scv
import matplotlib.pyplot as plt
import os
import numpy as np
import scipy.sparse
import utils 

# ================= 🔧 配置区域 =================
# 1.  .h5ad 数据集路径
H5AD_PATH = 'velocity/indicators/gene/DentateGyrus_lamanno_analyze/results_for_analyze/celldancer.h5ad' 

# 2. Branching Genes 列表
# gene_list=['Diaph3','Dpp6','Etv1','Nnat','Nr4a2','Tubb3','ANK3']
# gene_list=['Malat1','Dab1','Mmd2','Plppr4','Ppp2r2b','Rasgrp1','Runx1t1']
# gene_list=['Plpp3','Rbfox3']
gene_list=['Gpm6b','Ncald','Diaph3','Dab1','Dpp6','Malat1','Nnat','Tubb3','Plppr4','Mmd2',]#最终选定这10个

BRANCHING_GENES = gene_list

cluster_edges = [
    # # --- 核心主干：从干细胞出发 ---
    ('RadialGlia', 'nIPC'),       
    ('nIPC', 'GlialProg'),  

    # # --- 神经元分支 (Neuron Lineage) ---
    ('Nbl2', 'ImmGranule1'),
    ('ImmGranule1', 'ImmGranule2'),
    ('ImmGranule2', 'Granule'),
    ('Nbl2', 'CA'),
    ('CA', 'CA1-Sub'),
    ('CA', 'CA2-3-4'),

    # # --- 胶质细胞分支 (Glia Lineage) ---
    ('GlialProg', 'OPC'),
]

# 3. 颜色标签
COLOR_KEY = 'clusters' 
plt.rcParams['figure.facecolor'] = 'white'

# 4. 图片输出目录
OUTPUT_DIR = 'velocity/indicators/gene/DentateGyrus_lamanno_analyze/branching_genes/results'

# 5. [新增] 分数保存路径 (txt)
SCORE_FILE_PATH = 'velocity/indicators/gene/DentateGyrus_lamanno_analyze/branching_genes/score.txt'
# ==============================================

def plot_celldancer_genes():
    
    # 1. 准备环境
    # 创建图片目录
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
    
    # [新增] 创建分数文件的父目录 (如果不存在)
    score_dir = os.path.dirname(SCORE_FILE_PATH)
    if not os.path.exists(score_dir):
        os.makedirs(score_dir)
        print(f"已创建分数保存目录: {score_dir}")

    # 如果需要，可以在这里先清空旧的 score.txt，或者选择保留追加
    # with open(SCORE_FILE_PATH, 'w') as f: f.write("") # 取消注释此行可每次运行时清空旧文件

    scv.settings.set_figure_params('scvelo', dpi=100, frameon=False)

    # 2. 加载数据
    if not os.path.exists(H5AD_PATH):
        print(f"找不到数据集 {H5AD_PATH}，请修改脚本中的路径。")
        return

    # 提取方法名
    method_name = os.path.splitext(os.path.basename(H5AD_PATH))[0]
    print(f"检测到方法名称: {method_name}")

    print(f"正在加载数据集: {H5AD_PATH} ...")
    try:
        adata = sc.read_h5ad(H5AD_PATH)
        print(adata)
        if 'velocity_unspliced' in adata.layers:

            adata.layers['velocity_u']=adata.layers['velocity_unspliced']
        if 'du_dt' in adata.layers:

            adata.layers['velocity_u']=adata.layers['du_dt']
        if 'velocity_reconstructed' in adata.layers:

            adata.layers['velocity']=adata.layers['velocity_reconstructed']
        if 'velocity_gv' in adata.layers:

            adata.layers['velocity']=adata.layers['velocity_gv']
        if 'velo_u' in adata.layers:
            adata.layers[ 'velocity_u'] = adata.layers['velo_u' ] / np.sum(np.abs(adata.layers['velo_u' ]), 0)
        # # 🟢 细胞过滤模块
        # target_clusters = list(set([node for edge in cluster_edges for node in edge]))
        # print(f"\n正在执行细胞过滤，只保留以下类型: {target_clusters}")
        
        # if COLOR_KEY in adata.obs.columns:
        #     n_original = adata.n_obs
        #     adata = adata[adata.obs[COLOR_KEY].isin(target_clusters)].copy()
        #     n_filtered = adata.n_obs
        #     print(f"过滤完成: {n_original} -> {n_filtered} cells")
        #     if n_filtered == 0:
        #         print("❌ 错误: 过滤后没有剩余细胞！")
        #         return
        # else:
        #     print(f"❌ 警告: adata.obs 中找不到 '{COLOR_KEY}' 列。")

        # 稀疏矩阵处理
                #  1. 先把 NaN/inf 清掉
        vel = adata.layers['velocity'].astype(np.float64, copy=False)
        vel = np.nan_to_num(vel, nan=0.0, posinf=0.0, neginf=0.0)

    # 2. 
        vel = np.clip(vel, -1e+10, 1e+10)
        adata.layers['velocity'] = vel
                #  1. 先把 NaN/inf 清掉
        vel = adata.layers['velocity_u'].astype(np.float64, copy=False)
        vel = np.nan_to_num(vel, nan=0.0, posinf=0.0, neginf=0.0)

    # 2. 
        vel = np.clip(vel, -1e+10, 1e+10)
        adata.layers['velocity_u'] = vel
        layers_to_fix = ['spliced', 'unspliced', 'velocity', 'DeepKINET_velocity']
        print("正在检查并转换稀疏矩阵层...")
        for layer in layers_to_fix:
            if layer in adata.layers:
                if scipy.sparse.issparse(adata.layers[layer]):
                    print(f"  -> converting sparse layer: {layer}")
                    adata.layers[layer] = adata.layers[layer].toarray()
                else:
                    print(f"  -> layer '{layer}' is already dense.")

        if scipy.sparse.issparse(adata.X):
            adata.X = adata.X.toarray()
            
        scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
        scv.pp.neighbors(adata)
        
    except Exception as e:
        print(f"❌ 读取数据失败: {e}")
        return

    plot_color = COLOR_KEY if COLOR_KEY in adata.obs.columns else None

    # 3. 遍历基因并绘图、计算分数
    print(f"准备绘制 {len(BRANCHING_GENES)} 个基因...")
    
    success_count = 0
    results = []
    
    for gene in BRANCHING_GENES:
        if gene not in adata.var_names:
            print(f"  [跳过] 数据集中未找到基因: {gene}")
            continue
            
        print(f"  正在绘制: {gene} ...")
        fig,ax = plt.subplots(figsize=(10, 6))
        # 图片路径: 基因名_方法名.png
        save_path = os.path.join(OUTPUT_DIR, f"{gene}_{method_name}.png")
        print(adata.layers['velocity'][:, adata.var_names.get_loc(gene)].ravel())
        print(adata.layers['velocity_u'][:, adata.var_names.get_loc(gene)].ravel())
        utils.plot_phase_portrait(
            adata,
            gene=gene,           # 你想看的基因名
            ax=ax,
            
            # 指定数据层
            layer_s='Ms',           # X轴: Spliced Moments
            layer_u='Mu',           # Y轴: Unspliced Moments
            layer_v_s='velocity',   # dS/dt: 也就是 scVelo 计算的主速率
            layer_v_u='velocity_u',         # dU/dt: scVelo 动力学模式通常不输出这个，设为 None 即可
            
            # 视觉参数
            color_key='clusters',   # 用哪一列来标颜色 (在 adata.obs 中)
            arrow_grid=(25, 25)     # 箭头密度，数值越大箭头越密 (传入 sampling_neighbors)
        )
        plt.savefig(save_path, dpi=600, facecolor='white', transparent=False)

        plt.close(fig) # 关闭图片释放内存

        print(f"正在计算分数: {gene} ...")
        try:
            raw_scores, global_score = utils.cross_boundary_correctness_gene(
                adata, 
                gene_name=gene,
                k_cluster='clusters', 
                k_velocity='velocity', 
                cluster_edges=cluster_edges, 
                return_raw=False
            )
            
            print(f"  -> Gene {gene}: CBDir: {global_score:.4f}")
            
            # ==========================================
            # 🟢 新增：将分数写入 txt 文件
            # ==========================================
            # 格式： 基因名: GlobalScore | Raw: {Details}
            # 使用追加模式 'a'，防止覆盖之前的记录
            with open(SCORE_FILE_PATH, 'a', encoding='utf-8') as f:
                # 写入更详细的信息，既包含global也包含raw
                line_content = f"{method_name}:{gene}: {global_score:.5f} | Raw_Scores: {str(raw_scores)}\n"
                f.write(line_content)
                print(f"  -> 已保存分数至: {SCORE_FILE_PATH}")
            # ==========================================

            success_count += 1
            
        except Exception as e:
            print(f"  -> 计算出错: {e}")

    print(f"\n✅ 完成！成功处理 {success_count}/{len(BRANCHING_GENES)} 个基因。")
    print(f"📂 图片目录: '{OUTPUT_DIR}'")
    print(f"📝 分数文件: '{SCORE_FILE_PATH}'")

if __name__ == "__main__":
    plot_celldancer_genes()