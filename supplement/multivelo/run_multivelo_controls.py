"""Small paired RNA-only refits with inferred versus unrestricted phase direction."""
from pathlib import Path
import argparse,sys,json,hashlib,platform,inspect,time
import anndata as ad
import numpy as np
import pandas as pd
import multivelo as mv
import multivelo.dynamical_chrom_func as implementation
from scipy import sparse

p=argparse.ArgumentParser();p.add_argument('--repository-root',type=Path,required=True);p.add_argument('--artifact-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--all-genes',action='store_true');p.add_argument('--no-output-smoothing',action='store_true')
p.add_argument('--candidate-gate-control',action='store_true',help='Only remove the on-direction repression-candidate skip; retain native on penalties')
p.add_argument('--case-gene',type=str,help='Pilot one gene in the normal + branching realization')
p.add_argument('--branching-only',action='store_true',help='Both evaluated branching subsets, 20 genes each')
args=p.parse_args()
assert not args.output.exists();args.output.mkdir(parents=True)
sys.path.insert(0,str(args.repository_root))
from evaluation.gene.config import SIM_PATHS,SIM_CONFIGS
from evaluation.gene.run_reference30_gene_cbdir import score_gene,align_reference_neighbors

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def arr(x):return np.asarray(x.toarray() if sparse.issparse(x) else x,dtype=float)
def score(o,g,c,neighbors):
 j=o.var_names.get_loc(g);ms=arr(o.layers['Ms'][:,j]).ravel();mu=arr(o.layers['Mu'][:,j]).ravel();vs=arr(o.layers['velo_s'][:,j]).ravel();vu=arr(o.layers['velo_u'][:,j]).ravel()
 v=np.column_stack([vs,vu]);gt=np.column_stack([arr(o.layers['velocity_s_true'][:,j]).ravel(),arr(o.layers['velocity_u_true'][:,j]).ravel()]);nv=np.linalg.norm(v,axis=1);nt=np.linalg.norm(gt,axis=1);ok=(nv>1e-6)&(nt>1e-6)
 co=np.clip((v[ok]*gt[ok]).sum(1)/(nv[ok]*nt[ok]),-1,1).mean()
 cb=score_gene(ms,mu,vs,vu,o.obs['bin'].astype(str).to_numpy(),neighbors,c['edges'])[0]
 return float(cb),float(co)

if args.candidate_gate_control:
 import ast,copy,difflib,subprocess
 assert not args.no_output_smoothing
 assert bool(args.case_gene)!=bool(args.branching_only), 'Select a pilot gene or both branching subsets'
 native_distance=implementation.calculate_dist_and_time
 native_velocity=implementation.velocity
 native_fit_dyn=implementation.ChromatinDynamical.fit_dyn
 native_source=inspect.getsource(native_distance)
 target="elif direction == 'on':"
 assert native_source.count(target)==1
 modified_source=native_source.replace(target,"elif False:  # diagnostic: retain repression candidates for direction='on'")
 expected_ast=ast.parse(native_source)
 targets=[n for n in ast.walk(expected_ast) if isinstance(n,ast.If) and ast.dump(n.test)==ast.dump(ast.parse("direction == 'on'",mode='eval').body)]
 assert len(targets)==1
 targets[0].test=ast.Constant(value=False)
 assert ast.dump(expected_ast)==ast.dump(ast.parse(modified_source)), 'Unexpected source changes'
 namespace=implementation.__dict__.copy()
 exec(compile(modified_source,'<candidate-gate-only-diagnostic>','exec'),namespace)
 open_distance=namespace['calculate_dist_and_time']
 (args.output/'calculate_dist_and_time_candidate_gate.py').write_text(modified_source)
 (args.output/'candidate_gate_only.patch').write_text(''.join(difflib.unified_diff(native_source.splitlines(True),modified_source.splitlines(True),fromfile='installed_function',tofile='in_memory_function')))

 # Source-level change and small numerical checks precede any refitting.
 c=np.ones(8);u=np.array([.1,.5,1,2,3,1,.3,.1]);s=np.array([.1,.2,.5,1,2,3,4,5.])
 test_args=(c,u,s,np.array([.1,7.,20.]),.1,3.,1.,.5,1.,1.)
 checks={}
 for direction in ['off','complete']:
  left=native_distance(*test_args,model=2,direction=direction,rna_only=True,t=200)
  right=open_distance(*test_args,model=2,direction=direction,rna_only=True,t=200)
  for a,b in zip(left[:3],right[:3]):np.testing.assert_array_equal(a,b)
  checks[f'{direction}_path_unchanged']=True
 left=open_distance(*test_args,model=2,direction='on',rna_only=True,t=200,penalize_gap=False)
 right=native_distance(*test_args,model=2,direction='complete',rna_only=True,t=200,penalize_gap=False)
 for a,b in zip(left[:3],right[:3]):np.testing.assert_array_equal(a,b)
 checks['open_on_matches_both_candidates_without_complete_gap_penalty']=True
 checks['only_one_ast_condition_changed']=True
 (args.output/'patch_validation.json').write_text(json.dumps(checks,indent=2)+'\n')

 selection=({'branching':[args.case_gene]} if args.case_gene else
            {ds:[g for g in SIM_CONFIGS[ds]['genes'] if not g.startswith('OneSidedBoost')] for ds in ['branching','all']})
 parameters=dict(rna_only=True,max_iter=5,init_mode='invert',parallel=False,n_jobs=1,save_plot=False,n_anchors=500,extra_color_key=None,verbose=False)
 source_file=Path(inspect.getfile(implementation))
 expected_source='ac9a998a9b896ed142f0113b19416ef1e0afcf9cc1978487419ca842d9d34020'
 assert sha(source_file)==expected_source
 frozen={str(source_file):sha(source_file)}
 for ds in selection:
  for path in [args.artifact_root/SIM_PATHS[ds]['MultiVelo'],args.repository_root/f'evaluation/gene/results/reference_neighbors/{ds}.npz']:
   frozen[str(path)]=sha(path)
 for relative in ['evaluation/gene/config.py','evaluation/gene/run_reference30_gene_cbdir.py','evaluation/gene/results/gene_cbdir_reference30_all.csv']:
  path=args.repository_root/relative;frozen[str(path)]=sha(path)
 manifest=dict(repository_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=args.repository_root,text=True).strip(),
               multivelo=mv.__version__,python=platform.python_version(),numpy=np.__version__,parameters=parameters,
               selection=selection,input_hashes=frozen,script_sha256=sha(Path(__file__)),patch_validation=checks,
               intervention="Keep automatic direction='on' and all its native penalties; only disable repression-candidate exclusion inside calculate_dist_and_time",
               randomness='No added seed or randomized resampling; package GMM uses random_state=2021',
               baseline_gate='Default velocities/metrics must reproduce frozen outputs; initial parameter/graph fingerprints must match before comparing effects',
               primary_results_modified=False,installed_source_modified=False)
 (args.output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
 context={};initial={};final_state={}
 def fit_dyn_audit(self):
  key=(context['dataset'],context['condition'],self.gene)
  context['gene']=self.gene
  initial[key]=dict(direction=self.direction,partial=bool(self.partial),params=self.params.tolist(),
                    c_sha256=hashlib.sha256(np.ascontiguousarray(self.c).tobytes()).hexdigest(),
                    u_sha256=hashlib.sha256(np.ascontiguousarray(self.u).tobytes()).hexdigest(),
                    s_sha256=hashlib.sha256(np.ascontiguousarray(self.s).tobytes()).hexdigest())
  return native_fit_dyn(self)
 def velocity_audit(t,switch,state,*vargs,**vkwargs):
  if state is not None:
   final_state[(context['dataset'],context['condition'],context['gene'])]=(np.asarray(t).copy(),np.asarray(state).copy())
  return native_velocity(t,switch,state,*vargs,**vkwargs)
 implementation.ChromatinDynamical.fit_dyn=fit_dyn_audit
 implementation.velocity=velocity_audit
 rows=[];branch_rows=[];cell_rows=[];initial_rows=[];baseline_checks=[]
 try:
  for ds,genes in selection.items():
   path=args.artifact_root/SIM_PATHS[ds]['MultiVelo'];original=ad.read_h5ad(path)
   asset=np.load(args.repository_root/f'evaluation/gene/results/reference_neighbors/{ds}.npz')
   neighbors=align_reference_neighbors(original,asset['obs_names'].astype(str),asset['indices'])
   graph_fingerprint=None
   for condition,distance_function in [('refit_default',native_distance),('on_allow_repression',open_distance)]:
    implementation.calculate_dist_and_time=distance_function
    context.update(dataset=ds,condition=condition)
    start=time.time();print('CANDIDATE_GATE',ds,condition,genes,flush=True)
    result=mv.recover_dynamics_chrom(original.copy(),adata_atac=None,gene_list=genes,**parameters)
    graph=result.obsp['_RNA_conn'].tocsr()
    fingerprint=hashlib.sha256(graph.data.tobytes()+graph.indices.tobytes()+graph.indptr.tobytes()).hexdigest()
    if graph_fingerprint is None:graph_fingerprint=fingerprint
    else:assert fingerprint==graph_fingerprint,'Neighbor graph changed between conditions'
    labels=result.obs['bin'].astype(str).to_numpy();branches=np.array([x.rsplit('_',1)[0] for x in labels])
    for g in genes:
     key=(ds,condition,g);j=result.var_names.get_loc(g);oj=original.var_names.get_loc(g)
     init=initial[key]
     assert init['direction']=='on' and result.var.loc[g,'fit_direction']=='on'
     initial_rows.append(dict(dataset=ds,condition=condition,gene=g,**init))
     if condition=='on_allow_repression':
      assert init==initial[(ds,'refit_default',g)],'Initialization changed'
     cb,co=score(result,g,SIM_CONFIGS[ds],neighbors);bcb,bco=score(original,g,SIM_CONFIGS[ds],neighbors)
     velocities=np.column_stack([arr(result.layers[k][:,j]).ravel() for k in ['velo_s','velo_u']])
     archived=np.column_stack([arr(original.layers[k][:,oj]).ravel() for k in ['velo_s','velo_u']])
     if condition=='refit_default':
      error=float(np.max(np.abs(velocities-archived)))
      assert error<1e-9 and abs(cb-bcb)<1e-10 and abs(co-bco)<1e-10,(ds,g,error,cb-bcb,co-bco)
      baseline_checks.append(dict(dataset=ds,gene=g,max_abs_velocity_error=error,cbdir_error=cb-bcb,cossim_error=co-bco))
     truth=np.column_stack([arr(result.layers[k][:,j]).ravel() for k in ['velocity_s_true','velocity_u_true']])
     nv=np.linalg.norm(velocities,axis=1);nt=np.linalg.norm(truth,axis=1);valid=(nv>1e-6)&(nt>1e-6)
     cos=np.full(result.n_obs,np.nan);cos[valid]=np.clip((velocities[valid]*truth[valid]).sum(1)/(nv[valid]*nt[valid]),-1,1)
     native_time,native_state=final_state[key];off=native_state>=2
     fit=result.var.loc[g]
     rows.append(dict(dataset=ds,gene=g,family=g.rsplit('_',1)[0],condition=condition,direction=fit.fit_direction,
                      cbdir=cb,cossim=co,saved_cbdir=bcb,saved_cossim=bco,valid_cells=int(valid.sum()),
                      alpha=fit.fit_alpha,beta=fit.fit_beta,gamma=fit.fit_gamma,switch=fit.fit_t_sw2,
                      unspliced_scaling=fit.fit_rescale_u,fit_loss=fit.fit_loss,
                      final_repression_fraction=float(off.mean()),wall_seconds=time.time()-start))
     for branch in ['BranchA','BranchB']:
      mask=branches==branch;good=mask&valid
      branch_rows.append(dict(dataset=ds,gene=g,family=g.rsplit('_',1)[0],condition=condition,branch=branch,n_cells=int(mask.sum()),valid_cells=int(good.sum()),
                              cossim=float(cos[good].mean()),final_repression_fraction=float(off[mask].mean()),
                              opposite_direction_fraction=float(np.mean(cos[good]<0)),positive_true_v_u_fraction=float(np.mean(truth[mask,1]>1e-6))))
     if args.case_gene:
      for i,cell in enumerate(result.obs_names):
       cell_rows.append(dict(dataset=ds,gene=g,condition=condition,cell_id=cell,branch=branches[i],bin=labels[i],
                             spliced=float(result.layers['Ms'][i,j]),unspliced=float(result.layers['Mu'][i,j]),
                             final_velocity_time=float(native_time[i]),final_velocity_state=int(native_state[i]),
                             v_s=float(velocities[i,0]),v_u=float(velocities[i,1]),true_v_s=float(truth[i,0]),true_v_u=float(truth[i,1]),cossim=float(cos[i])))
    result.write_h5ad(args.output/f'{ds}_{condition}.h5ad',compression='gzip')
    pd.DataFrame(rows).to_csv(args.output/'paired_refit_scores.csv',index=False)
    pd.DataFrame(branch_rows).to_csv(args.output/'branch_scores.csv',index=False)
    (args.output/'initialization_audit.json').write_text(json.dumps(initial_rows,indent=2)+'\n')
    pd.DataFrame(baseline_checks).to_csv(args.output/'baseline_reproduction.csv',index=False)
    if cell_rows:pd.DataFrame(cell_rows).to_csv(args.output/'cell_velocity_states.csv',index=False)
    print('CANDIDATE_GATE_RESULT',ds,condition,pd.DataFrame(rows).query('dataset == @ds and condition == @condition')[['gene','cbdir','cossim','final_repression_fraction']].to_json(orient='records'),flush=True)
 finally:
  implementation.calculate_dist_and_time=native_distance
  implementation.velocity=native_velocity
  implementation.ChromatinDynamical.fit_dyn=native_fit_dyn
 assert all(sha(Path(path))==digest for path,digest in frozen.items()),'Frozen source or input changed'
 manifest.update(inputs_unchanged_after_refits=True,initialization_equal_between_conditions=True,neighbor_graph_equal_between_conditions=True,baseline_reproduction_passed=True,n_gene_pairs=sum(map(len,selection.values())),
                 outputs={f.name:sha(f) for f in args.output.iterdir() if f.is_file() and f.name!='manifest.json'})
 frame=pd.DataFrame(rows)
 frame.groupby(['dataset','condition'])[['cbdir','cossim','final_repression_fraction']].mean().reset_index().to_csv(args.output/'summary.csv',index=False)
 manifest['outputs']['summary.csv']=sha(args.output/'summary.csv')
 (args.output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
 print('CANDIDATE_GATE_DONE',frame.to_json(orient='records'),flush=True)
 sys.exit(0)

selection={'branching':['MultiForward_1','BoostALL_5'],'all':['MultiForward_1','BoostALL_5','OneSidedBoost_0'],'murk':['OneSidedBoost_0']}
if args.all_genes:
 selection={ds:config['genes'] for ds,config in SIM_CONFIGS.items()}
manifest=dict(multivelo=mv.__version__,python=platform.python_version(),numpy=np.__version__,parameters=dict(rna_only=True,max_iter=5,init_mode='invert',parallel=False,n_jobs=1,n_anchors=500),randomness='No global seed change; package partial-trajectory GMM fixes random_state=2021',sources={str(args.artifact_root/SIM_PATHS[ds]['MultiVelo']):sha(args.artifact_root/SIM_PATHS[ds]['MultiVelo']) for ds in selection},source_code=dict(path=inspect.getfile(implementation),sha256=sha(Path(inspect.getfile(implementation)))),selection=selection,script_sha256=sha(Path(__file__)))
conditions=[('refit_default',None),('refit_complete','complete')]
if args.no_output_smoothing:
 implementation.smooth_scale=lambda conn,vector:vector.copy()
 conditions=[('refit_no_output_smoothing',None)]
 manifest['intervention']='Replace module smooth_scale with identity in this process only; installed source files untouched'
(args.output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
rows=[]
for ds,genes in selection.items():
 path=args.artifact_root/SIM_PATHS[ds]['MultiVelo'];manifest['sources'][str(path)]=sha(path)
 original=ad.read_h5ad(path);asset=np.load(args.repository_root/f'evaluation/gene/results/reference_neighbors/{ds}.npz');neighbors=align_reference_neighbors(original,asset['obs_names'].astype(str),asset['indices'])
 for condition,direction in conditions:
  start=time.time();print('REFIT',ds,condition,genes,flush=True)
  result=mv.recover_dynamics_chrom(original.copy(),adata_atac=None,gene_list=genes,max_iter=5,init_mode='invert',parallel=False,n_jobs=1,save_plot=False,rna_only=True,fit=True,n_anchors=500,extra_color_key=None,verbose=False,direction=direction)
  print('Finished in',time.time()-start,flush=True)
  for g in genes:
   cb,co=score(result,g,SIM_CONFIGS[ds],neighbors);bcb,bco=score(original,g,SIM_CONFIGS[ds],neighbors);fit=result.var.loc[g]
   rows.append(dict(dataset=ds,gene=g,condition=condition,direction=fit.fit_direction,cbdir=cb,cossim=co,saved_cbdir=bcb,saved_cossim=bco,alpha=fit.fit_alpha,beta=fit.fit_beta,gamma=fit.fit_gamma,switch=fit.fit_t_sw2,unspliced_scaling=fit.fit_rescale_u,fit_loss=fit.fit_loss,wall_seconds=time.time()-start))
  pd.DataFrame(rows).to_csv(args.output/'paired_refit_scores.csv',index=False)
  result.var.to_csv(args.output/f'{ds}_{condition}_parameters.csv')
 assert sha(path)==manifest['sources'][str(path)]
manifest['inputs_unchanged_after_refits']=True
(args.output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(pd.DataFrame(rows).to_string(index=False))
