"""Saved-fit scVelo controls; no training and no source modification."""
from pathlib import Path
import argparse,json,hashlib,sys,platform
import numpy as np
import pandas as pd
import anndata as ad
import scvelo as scv
from scipy import sparse
from scvelo.tools._em_model_utils import get_divergence,get_vars,compute_divergence
from scvelo.core import SplicingDynamics

p=argparse.ArgumentParser();p.add_argument('--repository-root',type=Path,required=True);p.add_argument('--artifact-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
assert not a.output.exists();a.output.mkdir(parents=True)
sys.path.insert(0,str(a.repository_root))
from evaluation.gene.config import SIM_PATHS,SIM_CONFIGS
from evaluation.gene.run_reference30_gene_cbdir import score_gene,align_reference_neighbors
from evaluation.gene.simulation.reference import load_velocity_reference
def sha(path):
 return hashlib.sha256(path.read_bytes()).hexdigest()

def arr(x):
 return np.asarray(x.toarray() if sparse.issparse(x) else x,dtype=float)

def cosines(v,truth):
 nv=np.linalg.norm(v,axis=-1);nt=np.linalg.norm(truth,axis=-1);good=(nv>1e-6)&(nt>1e-6)
 result=np.zeros(nv.shape);result[good]=np.clip(np.sum(v[good]*truth[good],axis=-1)/(nv[good]*nt[good]),-1,1)
 return float(result[good].mean()),float(result[nt>1e-6].mean()),int(good.sum())
rows=[];checks=[];frozen={}
for ds,c in SIM_CONFIGS.items():
 print('CONTROL',ds,flush=True)
 path=a.artifact_root/SIM_PATHS[ds]['scVelo_dyn'];frozen[str(path)]=sha(path)
 full=ad.read_h5ad(path);obj=full[:,c['genes']].copy()
 conn=scv.utils.get_connectivities(obj) if hasattr(scv.utils,'get_connectivities') else None
 if conn is None:
  from scvelo.preprocessing.neighbors import get_connectivities
  conn=get_connectivities(obj)
 opt=dict(kernel_width=None,normalized=True,var_scale=True,reg_par=None,min_confidence=1e-2,constraint_time_increments=False,fit_steady_states=True,fit_basal_transcription=None,use_connectivities=True,time_connectivities=True,use_latent_time=None)
 opt.update(dict(obj.uns['recover_dynamics']))
 baseline=np.stack([arr(obj.layers['velocity']),arr(obj.layers['velocity_u'])],axis=-1)
 variants={'saved':baseline}
 for name,change in [('source_recomputed',{}),('no_steady_states',dict(fit_steady_states=False)),('no_time_smoothing',dict(time_connectivities=False)),('no_state_or_time_smoothing',dict(use_connectivities=False,time_connectivities=False)),('no_min_confidence',dict(min_confidence=0)),('no_steady_or_min_confidence',dict(fit_steady_states=False,min_confidence=0))]:
  settings={**opt,**change}
  try:
   captured={}
   def capture_state(frame,event,value):
    if frame.f_code is not compute_divergence.__code__:
     return None
    if event=='return' and frame.f_locals.get('mode')=='velocity':
     captured['on']=np.asarray(frame.f_locals['o']).copy()
     captured['off']=np.asarray(frame.f_locals['o_']).copy()
    return capture_state
   previous_trace=sys.gettrace()
   try:
    if name=='source_recomputed' and ds=='branching':sys.settrace(capture_state)
    vs,vu=get_divergence(obj,mode='velocity',**settings)
   finally:
    sys.settrace(previous_trace)
   if name=='source_recomputed' and ds=='branching':
    j=obj.var_names.get_loc('BoostALL_5')
    on=captured['on'][:,j];off=captured['off'][:,j]
    pd.DataFrame(dict(cell_id=obj.obs_names,gene='BoostALL_5',bin=obj.obs['bin'].astype(str).to_numpy(),final_state=np.where(on,'induction',np.where(off,'repression','steady_or_unassigned')))).to_csv(a.output/'phase_assignment_case.csv',index=False)
   variants[name]=np.stack([vs,vu],axis=-1)
  except IndexError as exc:
   checks.append(dict(dataset=ds,unsupported_control=name,error=str(exc)))
   print('Unsupported upstream combination:',name,str(exc),flush=True)
 checks.append(dict(dataset=ds,max_abs_saved_recomputed=float(np.max(abs(variants['source_recomputed']-baseline)))))
 alpha,beta,gamma,scale,t_switch=get_vars(obj)
 tau=arr(obj.layers['fit_tau'])
 for name,tt in [('force_induction_fit_tau',tau),('force_induction_fit_tau_smoothed',conn@tau)]:
  u,s=SplicingDynamics(alpha=alpha,beta=beta,gamma=gamma).get_solution(tt,stacked=False)
  vs=beta*u-gamma*s;vu=(alpha-beta*u)*scale
  # Keep the native lower clipping used in scVelo's mode='velocity'.
  vs=np.maximum(vs,-arr(obj.layers['Ms']));vu=np.maximum(vu,-arr(obj.layers['Mu']))
  variants[name]=np.stack([vs,vu],axis=-1)
 truth=np.stack(load_velocity_reference(ds,obj.obs_names,c['genes']),axis=-1)
 asset=np.load(a.repository_root/f'evaluation/gene/results/reference_neighbors/{ds}.npz');neighbors=align_reference_neighbors(obj,asset['obs_names'].astype(str),asset['indices'])
 labels=obj.obs['bin'].astype(str).to_numpy();ms,mu=arr(obj.layers['Ms']),arr(obj.layers['Mu'])
 for name,v in variants.items():
  for j,g in enumerate(c['genes']):
   co,fix,n=cosines(v[:,j],truth[:,j]);cb=score_gene(ms[:,j],mu[:,j],v[:,j,0],v[:,j,1],labels,neighbors,c['edges'])[0]
   rows.append(dict(dataset=ds,gene=g,gene_class='MURK' if g.startswith('OneSidedBoost') else 'branching',family=g.rsplit('_',1)[0],variant=name,cbdir=cb,cossim=co,fixed_cossim=fix,valid_cells=n))
 pd.DataFrame(rows).to_csv(a.output/'gene_controls.csv',index=False)
pd.DataFrame(checks).to_csv(a.output/'reproduction_checks.csv',index=False)
pd.DataFrame(rows).groupby(['dataset','gene_class','variant'],sort=False)[['cbdir','cossim','fixed_cossim']].mean().reset_index().to_csv(a.output/'summary.csv',index=False)
assert all(sha(Path(k))==v for k,v in frozen.items())
(a.output/'manifest.json').write_text(json.dumps(dict(input_hashes=frozen,input_hashes_unchanged_after_run=True,scvelo=scv.__version__,numpy=np.__version__,python=platform.python_version(),script_sha256=sha(Path(__file__)),randomness='None; no fitting; saved parameters only'),indent=2)+'\n')
print(pd.DataFrame(checks).to_string(index=False));print(pd.read_csv(a.output/'summary.csv').to_string(index=False))
