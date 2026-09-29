#!/usr/bin/env python3
"""Reconstruct saved branching slope rankings without replacing historical CSVs.

The numerical operations follow the preserved project scripts. DentateGyrus
groups and granule/pyramidal ranks are recovered from dated editor snapshots.
"""
import argparse
import hashlib
import json
import os
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
from scipy import stats

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--input', type=Path, required=True)
parser.add_argument('--input-id', required=True)
parser.add_argument('--dataset', choices=['retina', 'DentateGyrus_lamanno'], required=True)
parser.add_argument('--include-neuron-glia', action='store_true',
                    help='Also run the unfiltered Neuron/Glia task recovered from the December snapshot')
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
if args.include_neuron_glia and args.dataset != 'DentateGyrus_lamanno':
    parser.error('Additional branching comparisons require DentateGyrus_lamanno')
output = args.output.resolve()
output.mkdir(parents=True, exist_ok=False)
adata = sc.read(str(args.input))
if args.dataset == 'retina':
    adata.obs['clusters'] = adata.obs['Annotation']
raw = adata.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.05, min_disp=0.1)
hvg = adata.var_names[adata.var.highly_variable].tolist()
adata = raw.copy()
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
genes = list(set(hvg + adata.var_names.tolist()))
adata = raw[:, genes].copy()
scv.pp.filter_and_normalize(adata, min_shared_counts=20)
scv.pp.moments(adata, n_neighbors=30, n_pcs=30)
del raw

def frame(layer):
    matrix = adata.layers[layer]
    return pd.DataFrame(matrix.toarray() if hasattr(matrix, 'toarray') else matrix,
                        index=adata.obs_names, columns=adata.var_names)

ms, mu = frame('Ms'), frame('Mu')

def fit(mask, max_s, max_u):
    n = int(np.sum(mask))
    if n < 10:
        raise ValueError('A comparison population contains fewer than ten cells')
    x, y = ms.loc[mask].values / max_s, mu.loc[mask].values / max_u
    mx, my = np.mean(x, axis=0), np.mean(y, axis=0)
    xc, yc = x - mx, y - my
    denominator = np.sum(xc ** 2, axis=0)
    with np.errstate(divide='ignore', invalid='ignore'):
        slope = np.sum(xc * yc, axis=0) / denominator
        slope[denominator == 0] = 0.0
        intercept = my - slope * mx
        prediction = slope * x + intercept
        residual = np.sum((prediction - y) ** 2, axis=0)
        se = np.sqrt(residual / (n - 2) / denominator)
        se[denominator == 0] = np.inf
    error = stats.t.ppf(0.975, n - 2) * se
    return slope, slope - error, slope + error

def trend(mask, ranks):
    labels = adata.obs['clusters'][mask].map(ranks).fillna(0).values
    result = []
    for values in ms.loc[mask].values.T:
        rho = stats.spearmanr(labels, values)[0] if np.std(values) else 0.0
        result.append(0.0 if np.isnan(rho) else rho)
    return np.array(result)

if args.dataset == 'retina':
    tasks = [('PR_vs_RGC', ['Neuroblast', 'PR'], ['Neuroblast', 'RGC'],
              {'Neuroblast': 1, 'PR': 2}, {'Neuroblast': 1, 'RGC': 1})]
    rank_status = 'matches_retained_retina_source_including_constant_RGC_rank'
else:
    tasks = [
        ('Granule_vs_Pyramidal', ['Nbl2', 'ImmGranule1', 'ImmGranule2', 'Granule'],
         ['Nbl2', 'CA', 'CA1-Sub', 'CA2-3-4'],
         {'Nbl2': 1, 'ImmGranule1': 2, 'ImmGranule2': 3, 'Granule': 4},
         {'Nbl2': 1, 'CA': 2, 'CA1-Sub': 3, 'CA2-3-4': 3})]
    rank_status = 'groups_and_ranks_recovered_from_20260107_editor_snapshot'
    if args.include_neuron_glia:
        tasks.append(
            ('Neuron_vs_Glia',
             ['RadialGlia', 'nIPC', 'Nbl1', 'Nbl2', 'ImmGranule1', 'ImmGranule2',
              'Granule', 'CA', 'CA1-Sub', 'CA2-3-4'],
             ['RadialGlia', 'GlialProg', 'ImmAstro', 'OPC'], None, None))
task_metadata = []
for name, cloud_a, cloud_b, ranks_a, ranks_b in tasks:
    missing = (set(cloud_a) | set(cloud_b)) - set(adata.obs['clusters'])
    if missing:
        raise ValueError('Missing clusters: {}'.format(sorted(missing)))
    mask_a, mask_b = adata.obs['clusters'].isin(cloud_a), adata.obs['clusters'].isin(cloud_b)
    union = mask_a | mask_b
    max_s, max_u = ms.loc[union].max().values, mu.loc[union].max().values
    max_s[max_s == 0] = 1.0
    max_u[max_u == 0] = 1.0
    slope_a, min_a, max_a = fit(mask_a, max_s, max_u)
    slope_b, min_b, max_b = fit(mask_b, max_s, max_u)
    gap = np.maximum(min_a - max_b, 0) + np.maximum(min_b - max_a, 0)
    label_a, label_b = ('neuron', 'glia') if name == 'Neuron_vs_Glia' else ('granule', 'pyramidal')
    slope_column_a, slope_column_b = 'slope_' + label_a, 'slope_' + label_b
    rho_column_a, rho_column_b = 'rho_' + label_a, 'rho_' + label_b
    table = pd.DataFrame({'gene': adata.var_names, 'separation_gap': gap,
                          slope_column_a: slope_a, slope_column_b: slope_b})
    table.sort_values('separation_gap', ascending=False).to_csv(output / (name + '_slopes.csv'), index=False)
    if name == 'Neuron_vs_Glia':
        task_metadata.append({'task': name, 'cloud_a': cloud_a, 'cloud_b': cloud_b,
                              'configuration_status': 'groups_recovered_from_20251216_editor_snapshot',
                              'n_a': int(mask_a.sum()), 'n_b': int(mask_b.sum()),
                              'n_ranked': len(table), 'trend_filter': 'not_present_in_historical_task'})
        continue
    rho_a, rho_b = trend(mask_a, ranks_a), trend(mask_b, ranks_b)
    table = table.rename(columns={'separation_gap': 'slope_gap'})
    table['trend_validity'] = np.maximum(np.abs(rho_a), np.abs(rho_b))
    table[rho_column_a], table[rho_column_b] = rho_a, rho_b
    table = table[['gene', 'slope_gap', 'trend_validity', slope_column_a,
                   slope_column_b, rho_column_a, rho_column_b]]
    filtered = table[table.trend_validity > 0.1].sort_values('slope_gap', ascending=False)
    filtered.to_csv(output / (name + '_filtered.csv'), index=False)
    task_metadata.append({'task': name, 'cloud_a': cloud_a, 'cloud_b': cloud_b,
                          'configuration_status': rank_status,
                          'ranks_a': ranks_a, 'ranks_b': ranks_b,
                          'n_a': int(mask_a.sum()), 'n_b': int(mask_b.sum()),
                          'n_ranked': len(table), 'n_filtered': len(filtered)})
digest = hashlib.sha256()
with args.input.open('rb') as stream:
    for chunk in iter(lambda: stream.read(1024 * 1024), b''):
        digest.update(chunk)
import anndata, scipy
(output / 'manifest.json').write_text(json.dumps({
    'status': 'reconstruction_check_not_manuscript_reselection',
    'dataset': args.dataset, 'input_id': args.input_id, 'input_sha256': digest.hexdigest(),
    'rank_configuration_status': rank_status, 'tasks': task_metadata,
    'preprocessing': {'min_mean': 0.05, 'min_disp': 0.1, 'target_sum': 1e4,
                      'min_shared_counts': 20, 'n_top_genes': 2000,
                      'n_neighbors': 30, 'n_pcs': 30, 'gene_union_order': 'list(set(...))'},
    'slope_interval': {'t_quantile': 0.975, 'degrees_of_freedom': 'n-2'},
    'trend_threshold_strictly_greater_than': 0.1,
    'pythonhashseed': os.environ.get('PYTHONHASHSEED'), 'python': platform.python_version(),
    'versions': {m.__name__: m.__version__ for m in [np, pd, sc, scv, anndata, scipy]},
}, indent=2) + '\n')
