# Get the MURK genes for GO analysis Figure 3D and rank the change in slope

# import scvelo as scv
# import numpy as np
# import pandas as pd
# import scipy
# import scanpy as sc

# sc.settings.vector_friendly = False
# scv.set_figure_params( dpi=300, dpi_save = 300, frameon=False, figsize = (7,4), format='png',fontsize=25)

# #Read data


# # adata_tot = scv.read(dir_data + '/Gastrulation/erythroid_lineage.h5ad')
# # adata_tot = sc.read_h5ad(dir_data + "/Gastrulation/erythroid_lineage.h5ad")
# adata_tot = sc.read_h5ad("velocity/data/Gastrulation/erythroid_lineage.h5ad")

# # adata_tot.var_names_make_unique() 
# # Prepare object

# # Define object characteristics

# celltype_colours = [
#     "#f9decf",
#     "#c9a997",
#     "#C72228",
#     "#f77b59",
#     "#EF4E22"]

# population_names = ['Blood progenitors 1','Blood progenitors 2', 'Erythroid1', 'Erythroid2', 'Erythroid3']

# # Subset object

# adata = adata_tot.copy()
# adata = adata[(adata.obs['celltype'] == 'Erythroid1') | 
#                   (adata.obs['celltype'] == 'Erythroid2') |
#                   (adata.obs['celltype'] == 'Erythroid3') | 
#                   (adata.obs['celltype'] == 'Blood progenitors 1') |
#                   (adata.obs['celltype'] == 'Blood progenitors 2'),:].copy()

# adata.uns['celltype_colors'] = celltype_colours
# adata_temp = adata.copy()

# vecB1 = np.array(adata.obs.celltype == 'Blood progenitors 1')
# vecB2 = np.array(adata.obs.celltype == 'Blood progenitors 2')
# vecE1 = np.array(adata.obs.celltype == 'Erythroid1')
# vecE2 = np.array(adata.obs.celltype == 'Erythroid2')
# vecE3 = np.array(adata.obs.celltype == 'Erythroid3')

# population_vectors = [vecB1, vecB2, vecE1, vecE2, vecE3]

# # Preprocessing 

# sc.pp.normalize_total(adata, target_sum=1e4)
# sc.pp.log1p(adata)
# sc.pp.highly_variable_genes(adata, min_mean=0.05,min_disp=0.1)
# list_hvg = adata.var_names[adata.var.highly_variable]

# adata = adata_temp.copy()
# scv.pp.filter_and_normalize(adata, min_shared_counts = 20, n_top_genes = 2000)

# list_scv = adata.var_names

# list_out = [gene for gene in list_scv if gene not in list_hvg]

# list_tot = list(list_hvg)

# list_tot.extend(list_out)

# # Imputation

# adata = adata_temp.copy()
# adata = adata[:,list_tot]
# # adata.write('adata_for_velocity_mouse.h5')

# scv.pp.filter_and_normalize(adata)
# scv.pp.moments(adata)

# # Murk genes research

# # Linear regression model: compute a slope for each gene and each population in phase plot space and store them in two dataframes

# # initialize matrix for spliced (mat_s) and unspliced (unspliced) counts
# mat_s = pd.DataFrame(adata.layers['Ms'].copy())
# mat_u = pd.DataFrame(adata.layers['Mu'].copy())
# mat_s.index = adata.obs_names
# mat_s.columns = adata.var_names
# mat_u.index = adata.obs_names
# mat_u.columns = adata.var_names

# # initialize dataframes for slopes (df_sl) and slope error (df_ds)
# df_sl = pd.DataFrame(index = adata.var_names)
# df_ds = pd.DataFrame(index = adata.var_names)

# # compute slopes and errors
# for pop in population_names:

#     vec = adata.obs['celltype'] == pop
    
#     x = mat_s.loc[vec]
#     y = mat_u.loc[vec]
    
#     mx = np.mean(x)
#     my = np.mean(y)
#     x_mx = x - mx
#     y_my = y - my
#     xy = np.sum(x_mx * y_my, axis = 0)
#     X2 = np.sum(x_mx**2,axis = 0)
#     sl = xy / X2
#     inter = my - sl * mx
#     n = np.sum(vec)
                
#     df_sl[pop] = pd.DataFrame(sl)     
    
#     pred = sl * x + inter
#     sse = (pred - y)**2
#     SSE = np.sum(sse, axis = 0)
                
#     ds = np.sqrt(SSE/(n-2)/X2)
#     df_ds[pop] = pd.DataFrame(ds)

# # initialize dataframes for slopes (df_exp), expression error (df_de) and 95% quantile (df_quant)
# df_exp = pd.DataFrame(index = adata.var_names)
# df_de = pd.DataFrame(index = adata.var_names)
# df_quan = pd.DataFrame(index = adata.var_names)

# # compute average expression and its error
# for pop in population_names:

#     vec = adata.obs['celltype'] == pop
#     n = np.sum(vec)
#     expr = np.mean(mat_s.loc[vec], axis = 0)
#     dex = np.std(mat_s.loc[vec], axis = 0)/np.sqrt(n)
    
#     df_exp[pop] = pd.DataFrame(expr)
#     df_de[pop] = pd.DataFrame(dex)

# for pop in population_names:
#     vec = adata.obs['celltype'] == pop
#     n = np.sum(vec)
    
#     df_quan[pop] = scipy.stats.t.ppf(0.95, n - 1)

# # compute minimum and maximum estimates of the slopes
# df_conf = df_ds * df_quan

# df_min = df_sl - df_conf
# df_max = df_sl + df_conf

# # decide if a change is slope is significative (Erythroid 3 expressed more than the others, its minimum slope grater than the previous maximum slope and positive)

# vec_mean = (((df_exp['Erythroid3'] >  df_exp['Erythroid2'])
#             |(df_exp['Erythroid3'] >  df_exp['Erythroid1']))
#             &(df_exp['Erythroid3'] >  df_exp['Blood progenitors 1'])
#             &(df_exp['Erythroid3'] >  df_exp['Blood progenitors 2']))

# vec_sl = (df_sl['Erythroid3'] > df_sl['Erythroid2']) 

# vec_pos =  df_sl['Erythroid3']>0

# vec_test = df_min['Erythroid3'] > df_max['Erythroid2']

# vec_tot = ((np.array(vec_sl) & np.array(vec_mean) & np.array(vec_pos) & np.array(vec_test)) | 
#  (vec_mean & np.array(df_sl['Erythroid3']<0)) )

# # How many MURK genes
# np.sum(vec_tot)

# # MURK genes for GO analysis Figure 3D
# murk_genes = adata.var_names[vec_tot]
# pd.DataFrame(murk_genes).to_csv('murk_genes_mouse.csv', index=None,header=None)

# # Recalculate the slopes scaling for the avarage of gene expression in order to range the genes for "MURKiness"

# df_sl2 = pd.DataFrame(index=adata.var_names)

# for pop in population_names:
#     print(pop)
#     vec = adata.obs['celltype'] == pop
    
#     x = mat_s.loc[vec]/np.max(mat_s.loc[vec])
#     y = mat_u.loc[vec]/np.max(mat_u.loc[vec])
    
#     mx = np.mean(x)
#     my = np.mean(y)
#     x_mx = x - mx
#     y_my = y - my
#     xy = np.abs(np.sum(x_mx * y_my, axis = 0))
#     X2 = np.sum(x_mx**2,axis = 0)
#     sl = xy / X2

#     df_sl2[pop] = pd.DataFrame(sl)

# # Rank genes: higher difference in scaled slope ranks higher
# df_rank = pd.DataFrame((df_min['Erythroid3'] - df_max['Erythroid2']))

# vec1 = np.array((df_max['Erythroid3']<0) & (df_min['Erythroid2']<0))

# df_rank.loc[vec1] = pd.DataFrame(-df_max.loc[vec1]['Erythroid3'] + df_min.loc[vec1]['Erythroid2'])

# vec2 = np.array((df_max['Erythroid3']<0) & (df_min['Erythroid2']>0))

# df_rank.loc[vec2] = pd.DataFrame(-df_max.loc[vec2]['Erythroid3'] + df_min.loc[vec2]['Erythroid2'])


# df_rank.columns = ['\u0394'+'m']
# df_rank['scaled \u0394'+'m'] = np.abs(df_sl2['Erythroid3'] - df_sl2['Erythroid2'])



# df_fin = df_rank.loc[murk_genes].sort_values(by = 'scaled \u0394'+'m',ascending = False)

# save_path = "saved/"
# df_fin.to_csv(save_path + 'test.csv') # Save final ranked MURK genes
