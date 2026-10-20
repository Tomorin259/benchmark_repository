#!/usr/bin/env python3
"""Compare current Figure 6 panels with an external archived panel directory."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--old', type=Path, required=True)
    parser.add_argument('--new', type=Path, required=True)
    parser.add_argument('--legacy-artifact-root', type=Path, help='Validate preserved renderer trajectory inputs separately from stored plot_data')
    parser.add_argument('--metrics-root', type=Path, help='Current evaluator output; required with --legacy-artifact-root')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--legacy-config', type=Path, help='Archived quantitative config, required for legacy-source validation')
    args = parser.parse_args()
    if args.legacy_artifact_root and (args.metrics_root is None or args.legacy_config is None):
        parser.error('--metrics-root and --legacy-config are required with --legacy-artifact-root')
    out = args.output
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(f'Output directory must be empty: {out}')
    out.mkdir(parents=True, exist_ok=True)
    olddata = args.old / 'plot_data'
    newdata = args.new / 'plot_data'
    summaries = []
    deltas = []
    means = []
    for newfile in sorted(newdata.glob('*.csv')):
        oldfile = olddata / newfile.name
        a, b = pd.read_csv(oldfile), pd.read_csv(newfile)
        if 'heatmap' in newfile.name:
            a = a.melt(id_vars='Method', var_name='Gene', value_name='Value')
            b = b.melt(id_vars='Method', var_name='Gene', value_name='Value')
            keys, values = ['Method', 'Gene'], ['Value']
        elif 'trajectories' in newfile.name:
            keys, values = ['Gene', 'Source', 'Relative_Position'], ['Display_Z']
        elif 'gene_metrics' in newfile.name:
            keys, values = ['Gene', 'Method', 'Parameter'], ['Correlation']
            for (method, parameter), rows in b.groupby(['Method', 'Parameter']):
                oldrows = a[(a.Method == method) & (a.Parameter == parameter)]
                oldmean, newmean = oldrows.Correlation.mean(), rows.Correlation.mean()
                means.append(dict(File=newfile.name, Method=method, Parameter=parameter,
                    Old_N=len(oldrows), New_N=len(rows), Old_Mean=oldmean,
                    New_Mean=newmean, Delta=newmean-oldmean,
                    Old_Display=f'{oldmean:.2f}', New_Display=f'{newmean:.2f}'))
        elif 'membership' in newfile.name:
            keys = ['Class', 'Gene']; values = [c for c in b.columns if c not in keys]
        else:
            keys = ['Class']; values = [c for c in b.columns if c not in keys]
        merged = a.merge(b, on=keys, how='outer', suffixes=('_old', '_new'), indicator=True, validate='one_to_one')
        shared = merged[merged._merge == 'both']
        maxdiff, changed, displayed = 0., 0, 0
        for column in values:
            x, y = pd.to_numeric(shared[column+'_old']).astype(float), pd.to_numeric(shared[column+'_new']).astype(float)
            delta = y-x
            finite = x.notna() & y.notna()
            mask = ((delta.abs() > 1e-12) & finite) | (x.isna() != y.isna())
            displaymask = (x.round(2) != y.round(2)) & ~(x.isna() & y.isna())
            changed += int(mask.sum()); displayed += int(displaymask.sum())
            if finite.any():
                maxdiff = max(maxdiff, float(delta[finite].abs().max()))
            for ix in shared.index[mask]:
                row = {k: shared.loc[ix,k] for k in keys}
                row.update(File=newfile.name, Value_Column=column, Old=shared.loc[ix,column+'_old'], New=shared.loc[ix,column+'_new'], Delta=delta.loc[ix])
                deltas.append(row)
        added = merged[merged._merge == 'right_only'][keys]
        removed = merged[merged._merge == 'left_only'][keys]
        if len(added): added.to_csv(out/(newfile.stem+'_added.csv'), index=False)
        if len(removed): removed.to_csv(out/(newfile.stem+'_removed.csv'), index=False)
        summaries.append(dict(File=newfile.name, Old_Rows=len(a), New_Rows=len(b), Shared_Rows=len(shared), Added_Rows=len(added), Removed_Rows=len(removed), Changed_Values_Above_1e_minus12=changed, Changed_At_2dp=displayed, Max_Absolute_Delta=maxdiff))
    pd.DataFrame(summaries).to_csv(out/'numeric_summary.csv', index=False)
    pd.DataFrame(deltas).to_csv(out/'numeric_differences.csv', index=False)
    pd.DataFrame(means).to_csv(out/'organoid_selected_gene_means.csv', index=False)
    if args.legacy_artifact_root:
        cfg = json.loads(args.legacy_config.read_text())
        sources = []
        display_names = {'Dynamo': 'Dynamo(label data)'}
        for parameter in ('alpha', 'gamma'):
            oldplot = pd.read_csv(olddata / f'rpe1_{parameter}_trajectories.csv')
            newplot = pd.read_csv(newdata / f'rpe1_{parameter}_trajectories.csv')
            for method, relative in cfg['inputs']['rpe1']['trajectories'][parameter].items():
                oldmatrix = pd.read_csv(args.legacy_artifact_root / relative, index_col=0)
                newmatrix = pd.read_csv(args.metrics_root / 'trajectories' / 'rpe1' / f'{method}_{parameter}.csv', index_col=0)
                newmatrix = newmatrix.reindex(oldmatrix.index)
                matrix_delta = np.nanmax(np.abs(oldmatrix.to_numpy(float)-newmatrix.to_numpy(float)))
                old_max, new_max = 0., 0.
                for gene in cfg['trajectory_genes']:
                    if gene not in oldmatrix.index:
                        continue
                    v = pd.Series(oldmatrix.loc[gene].to_numpy(float))
                    z = (v-v.mean())/v.std() if v.std() else v-v.mean()
                    for records, target in [(oldplot, 'old'), (newplot, 'new')]:
                        values = records[(records.Gene == gene) & (records.Source == display_names.get(method,method))].Display_Z.to_numpy(float)
                        delta = float(np.nanmax(np.abs(z.to_numpy()-values)))
                        if target == 'old': old_max=max(old_max,delta)
                        else: new_max=max(new_max,delta)
                sources.append(dict(Parameter=parameter, Method=method, Preserved_Static_Path=relative,
                    Max_Static_vs_Central_Processed_Delta=matrix_delta,
                    Max_Static_Z_vs_Preserved_PlotData_Delta=old_max,
                    Max_Static_Z_vs_Current_PlotData_Delta=new_max))
        pd.DataFrame(sources).to_csv(out/'trajectory_source_validation.csv', index=False)
    groups = {
        'A': [(f'A {method}', f'figure5A_rpe1/{stem}.png', f'figure6A_rpe1/{stem}.png') for method,stem in [('cellDancer','celldancer'),('DeepVelo','deepvelo'),('DeepKINET','deepkinet'),('TSvelo','tsvelo')]],
        'E': [(f'E {method}', f'figure5E_organoid/{stem}_organoid_stream_title_removed.png', f'figure6E_organoid/{stem}_organoid_stream_title_removed.png') for method,stem in [('cellDancer','celldancer'),('DeepKINET','DeepKINET'),('DeepVelo','deepvelo')]],
        'quantitative': [(label, stem+'.png', stem+'.png') for label,stem in [('B alpha','rpe1_alpha_trajectories'),('B gamma','rpe1_gamma_trajectories'),('C','rpe1_alpha_gamma_heatmap'),('D cellDancer','rpe1_celldancer_kinetic_classes'),('D DeepVelo','rpe1_deepvelo_kinetic_classes'),('F branch 1','organoid_branch1_boxplot'),('F branch 2','organoid_branch2_boxplot')]]}
    panelrecords=[]
    for group, panels in groups.items():
        fig, axes = plt.subplots(len(panels), 2, figsize=(16, 4*len(panels)), squeeze=False)
        for row,(label,oldrel,newrel) in enumerate(panels):
            paths=[args.old/oldrel,args.new/newrel]
            for column,path in enumerate(paths):
                axes[row,column].imshow(mpimg.imread(path)); axes[row,column].axis('off')
                axes[row,column].set_title(label+(' — archived panel' if column==0 else ' — current central result'))
            panelrecords.append(dict(Panel=label, Old=oldrel, New=newrel, Old_SHA256=sha256(paths[0]), New_SHA256=sha256(paths[1])))
        fig.tight_layout();fig.savefig(out/(group+'_old_new.jpg'),dpi=150);plt.close(fig)
    pd.DataFrame(panelrecords).to_csv(out/'panel_hashes.csv',index=False)
    lines=['# Figure 6 panel comparison','','Numerical tables and panel images are compared between the supplied --old and --new directories. numeric_differences.csv records changed values; panel_hashes.csv identifies the compared images.','','## Numeric comparison','', '| Data | Old / new rows | Added / removed | Changed at 2 decimals | Max absolute delta |','| --- | --- | --- | --- | --- |']
    for row in summaries:
        lines.append(f"| {row['File']} | {row['Old_Rows']} / {row['New_Rows']} | {row['Added_Rows']} / {row['Removed_Rows']} | {row['Changed_At_2dp']} | {row['Max_Absolute_Delta']:.9g} |")
    if means:
        lines += ['', '[Selected-gene mean comparison](organoid_selected_gene_means.csv)']
    if args.legacy_artifact_root:
        lines += ['', '[Trajectory source comparison](trajectory_source_validation.csv)']
    lines += ['', '## Images', '', '- [RPE1 streams A](A_old_new.jpg)', '- [Organoid streams E](E_old_new.jpg)', '- [Quantitative panels B–D and F](quantitative_old_new.jpg)']
    (out/'README.md').write_text('\n'.join(lines)+'\n')
    print(pd.DataFrame(summaries).to_string(index=False))
    print(pd.DataFrame(means).to_string(index=False))

if __name__ == '__main__': main()
