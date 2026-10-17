"""Inspect the exact slope/thickness decision inputs, without fitting or edits."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import anndata as ad

repo=Path(sys.argv[1]);artifact=Path(sys.argv[2]);out=Path(sys.argv[3])
sys.path.insert(0,str(repo))
from evaluation.gene.config import SIM_CONFIGS,SIM_PATHS
rows=[]
for ds,cfg in SIM_CONFIGS.items():
 a=ad.read_h5ad(artifact/SIM_PATHS[ds]['MultiVelo'])
 for gene in cfg['genes']:
  j=a.var_names.get_loc(gene)
  u=np.asarray(a.layers['Mu'][:,j]).ravel().astype(float);s=np.asarray(a.layers['Ms'][:,j]).ravel().astype(float)
  u-=u.min();s-=s.min()
  keep=((u>0)|(s>0))&(u<=np.percentile(u,99.8))&(s<=np.percentile(s,99.8))
  u,s=u[keep],s[keep]
  active=(u>=0.1*u.max())&(s>=0.1*s.max());ua,sa=u[active],s[active]
  high=(u>=np.percentile(ua,95))|(s>=np.percentile(sa,95))
  slope=float(u[high]@s[high]/(s[high]@s[high]))
  delta=ua-slope*sa
  on=float(np.sum(delta[delta>=0]**2));off=float(np.sum(delta[delta<0]**2))
  thickness=float(np.linalg.svd(np.column_stack([ua/u.max(),sa/s.max()]),compute_uv=False)[1])
  gate='on_SSE_gt_10x_off' if on>10*off else 'off_SSE_gt_10x_on' if off>10*on else 'depends_on_partial_and_thickness'
  rows.append(dict(dataset=ds,gene=gene,family=gene.rsplit('_',1)[0],initial_slope=slope,on_squared_residual=on,off_squared_residual=off,on_off_ratio=on/off if off>0 else np.inf,thickness=thickness,narrow=thickness<1.5,direct_gate=gate,saved_direction=a.var.loc[gene,'fit_direction']))
x=pd.DataFrame(rows);x.to_csv(out,index=False)
print(x.groupby(['dataset','family','direct_gate','saved_direction']).size())
print(x.groupby(['dataset','family'])[['initial_slope','on_off_ratio','thickness']].agg(['min','max']).round(4).to_string())
