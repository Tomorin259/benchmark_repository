#!/usr/bin/env python3
"""Export static single-gene displays from frozen inputs, without training."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from export_data import export_data, GENE_IDS, SIM_IDS, KINETICS


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def display_values(values, zscore=False):
    """Figure 6 RPE1 display: sample SD, finite values only, no gap filling."""
    import numpy as np
    import pandas as pd
    series = pd.Series(values, dtype=float)
    if zscore:
        sd = series.std()
        series = series - series.mean()
        if pd.notna(sd) and sd != 0:
            series /= sd
    return [float(v) if np.isfinite(v) else None for v in series]


def tsvelo_organoid(artifact_root, output):
    """Materialize the same branch processing and verify every saved score."""
    import anndata as ad
    import numpy as np
    import pandas as pd
    from evaluation.kinetics_recovery.expanded_reference import frozen_inputs, score_tsvelo
    from evaluation.kinetics_recovery.metrics import (
        calculate_metrics, csv_roundtrip, layer_frame, organoid_raw_from_layer,
        process_organoid, read_organoid_reference, select_organoid_branch)
    module = ROOT / 'evaluation/kinetics_recovery'
    inputs = frozen_inputs(artifact_root)
    # This existing evaluator checks source hashes, gene IDs and overlapping cells.
    scores, mapping = score_tsvelo(inputs)
    committed = pd.read_csv(module / 'results/expanded_reference/tsvelo_organoid_gene_metrics.csv')
    committed['Undefined_reason'] = committed.Undefined_reason.fillna('')
    keys = ['Method', 'Dataset', 'Branch', 'Parameter', 'Gene']
    pd.testing.assert_frame_equal(scores.sort_values(keys).reset_index(drop=True),
                                  committed.sort_values(keys).reset_index(drop=True),
                                  check_dtype=False, atol=1e-12, rtol=1e-12)
    settings = json.loads((module / 'config.json').read_text())['datasets']['organoid']
    rename = dict(zip(mapping.model_symbol, mapping.reference_symbol))
    obs_keys = [settings[k] for k in ('branch_key', 'cluster_key', 'pseudotime_key')]
    frames, metadata = [], []
    for key in ('tsvelo_lineage0', 'tsvelo_lineage1'):
        data = ad.read_h5ad(inputs[key])
        frames.append(layer_frame(data[:, list(rename)], 'alpha').rename(columns=rename))
        metadata.append(data.obs[obs_keys].copy())
    alpha = pd.concat(frames).groupby(level=0).mean()
    meta = pd.concat(metadata)
    meta = meta[~meta.index.duplicated(keep='first')].loc[alpha.index].copy()
    for key in obs_keys[:2]:
        meta[key] = meta[key].astype(str).str.replace(r'\.0$', '', regex=True)
    merged = ad.AnnData(obs=meta, var=pd.DataFrame(index=alpha.columns))
    merged.layers['alpha'] = alpha.to_numpy()
    traces = {}
    for branch, branch_settings in settings['branches'].items():
        selected, time = select_organoid_branch(merged, settings, branch)
        raw = csv_roundtrip(organoid_raw_from_layer(selected, time, 'alpha', settings))
        prediction = csv_roundtrip(process_organoid(raw, settings))
        truth = csv_roundtrip(read_organoid_reference(inputs['table_s4'], branch, 'alpha',
                              branch_settings['reference_reverse_time']).loc[prediction.index])
        actual = calculate_metrics(prediction, truth).set_index('Gene')
        expected = committed[committed.Branch.eq(int(branch))].set_index('Gene')
        for column in ('Spearman_Corr', 'Cosine_Similarity', 'Valid_Samples'):
            np.testing.assert_allclose(actual.loc[expected.index, column], expected[column],
                                       atol=1e-12, rtol=1e-12, equal_nan=True)
        path = output / f'TSvelo_alpha_branch{branch}.csv'
        prediction.to_csv(path)
        traces[str(branch)] = (prediction, path)
    return traces


def export_visuals(output, gene_root=None, kinetics_root=None):
    import numpy as np
    import pandas as pd
    output = Path(output).resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError('Use an empty output directory; previous displays must remain recoverable')
    if output == ROOT or output in ROOT.parents or (output.is_relative_to(ROOT) and
                                                   not output.is_relative_to(ROOT / 'site')):
        raise ValueError('Repository outputs must stay under site/')
    output.mkdir(parents=True, exist_ok=True)
    scores = export_data()
    sources, files = {}, []

    def source(path):
        relative = str(path.relative_to(ROOT))
        if relative not in sources:
            sources[relative] = {'path': relative, 'sha256': sha256(path)}
        return sources[relative]

    def register(path):
        files.append({'path': str(path.relative_to(output)), 'sha256': sha256(path),
                      'bytes': path.stat().st_size})

    central = ROOT / KINETICS
    frozen = json.loads((central / 'manifest.json').read_text())['outputs_sha256']
    source(central / 'manifest.json')
    source(ROOT / 'evaluation/kinetics_recovery/plotting/figure6.py')
    source(ROOT / 'site/export_visuals.py')
    for relative, digest in frozen.items():
        if '/trajectories/' in relative or '/reference_' in relative:
            if sha256(ROOT / relative) != digest:
                raise ValueError(f'Frozen matrix hash mismatch: {relative}')
    extra = {}
    if kinetics_root:
        trace_dir = output / 'processed'
        trace_dir.mkdir()
        extra = tsvelo_organoid(Path(kinetics_root), trace_dir)
        for _, path in extra.values():
            register(path)
        source(ROOT / 'evaluation/kinetics_recovery/results/expanded_reference/manifest.json')
        source(ROOT / 'evaluation/kinetics_recovery/results/expanded_reference/tsvelo_organoid_gene_metrics.csv')

    curves = []
    ids = {'cellDancer': 'C15', 'DeepVelo': 'C16', 'DeepKINET': 'C13',
           'TSvelo': 'C14', 'Dynamo': 'dynamo_labelled'}
    for dataset, branch in [('rpe1', ''), ('organoid', '1'), ('organoid', '2')]:
        folder = Path(dataset) / (f'branch{branch}' if branch else '')
        for parameter in ('alpha', 'gamma'):
            reference_path = central / 'matrices' / folder / f'reference_{parameter}.csv'
            reference = pd.read_csv(reference_path, index_col=0)
            matrices = [('reference', reference, source(reference_path))]
            for path in sorted((central / 'trajectories' / folder).glob(f'*_{parameter}.csv')):
                matrices.append((ids[path.stem.rsplit('_', 1)[0]], pd.read_csv(path, index_col=0), source(path)))
            if dataset == 'organoid' and parameter == 'alpha' and branch in extra:
                frame, path = extra[branch]
                matrices.append(('C14', frame, {'path': f'visuals/{path.relative_to(output)}', 'sha256': sha256(path)}))
            eligible = {r['gene'] for r in scores['kinetics'] if r['dataset'] == f'sceu_{dataset}'
                        and r['branch'] == branch and r['parameter'] == parameter}
            for gene in reference.index:
                if gene not in eligible:
                    continue
                series = []
                for method, matrix, provenance in matrices:
                    if gene not in matrix.index:
                        continue
                    series.append({'method': method, 'x': np.linspace(0, 1, matrix.shape[1]).tolist(),
                                   'positions': list(matrix.columns),
                                   'scale': 'display_z' if dataset == 'rpe1' else 'reference' if method == 'reference' else 'prediction',
                                   'y': display_values(matrix.loc[gene], dataset == 'rpe1'), 'source': provenance})
                path = output / 'curves' / folder / parameter / f'{gene}.json'
                path.parent.mkdir(parents=True, exist_ok=True)
                record = {'dataset': f'sceu_{dataset}', 'branch': branch, 'parameter': parameter,
                          'gene': gene, 'transform': 'Display z-score (sample SD)' if dataset == 'rpe1'
                          else 'Saved profiles: prediction baseline ratio; reference Table S4 values', 'series': series}
                path.write_text(json.dumps(record, allow_nan=False, separators=(',', ':')))
                register(path)
                curves.append({k: record[k] for k in ('dataset', 'branch', 'parameter', 'gene')} |
                              {'path': f'visuals/{path.relative_to(output)}'})
    print(f'Exported {len(curves)} single-gene trajectories', flush=True)

    phases = []
    if gene_root:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from evaluation.gene.plot_gene_assets import load_scores, read_phase_inputs, ORDER, CASES
        from evaluation.gene.plotting import plot_phase_portrait
        from evaluation.gene.config import SIM_CONFIGS
        table = load_scores()
        source(ROOT / 'evaluation/gene/plot_gene_assets.py')
        source(ROOT / 'evaluation/gene/plotting/phase_portrait.py')
        source(ROOT / 'evaluation/gene/results/gene_cbdir_reference30_all.csv')
        for case in CASES:
            genes = sorted(table.loc[table.dataset.eq(case), 'gene'].unique())
            adatas, inputs = read_phase_inputs(case, table, Path(gene_root), genes)
            for gene in genes:
                dataset = SIM_IDS.get(case, case.lower())
                for method, data, provenance in zip(ORDER, adatas, inputs):
                    entry = {'dataset': dataset, 'gene': gene, 'method': GENE_IDS[method],
                             'source': provenance, 'path': None, 'status': 'gene_absent'}
                    if gene in data.var_names:
                        fig, ax = plt.subplots(figsize=(4, 4))
                        _, diagnostics = plot_phase_portrait(data, gene, ax=ax, legend_loc=None,
                             color_key='bin' if case in SIM_CONFIGS else 'clusters', return_diagnostics=True)
                        ax.set_title('')
                        ax.set_box_aspect(1)
                        fig.tight_layout()
                        path = output / 'phases' / dataset / gene / f'{GENE_IDS[method]}.png'
                        path.parent.mkdir(parents=True, exist_ok=True)
                        fig.savefig(path, dpi=150, facecolor='white')
                        plt.close(fig)
                        register(path)
                        entry.update(path=f'visuals/{path.relative_to(output)}', status='available', diagnostics=diagnostics)
                    phases.append(entry)
            print(f'Exported phase portraits: {case}, {len(genes)} genes', flush=True)
            del adatas
    for item in sources.values():
        if sha256(ROOT / item['path']) != item['sha256']:
            raise ValueError(f"Source changed during export: {item['path']}")
    index = {'schema_version': 1, 'source_commit': subprocess.check_output(
             ['git', '-C', str(ROOT), 'rev-parse', 'HEAD'], text=True).strip(),
             'curves': curves, 'phases': phases, 'sources': list(sources.values()), 'files': files}
    (output / 'index.json').write_text(json.dumps(index, allow_nan=False, separators=(',', ':')))
    return index


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--gene-artifact-root', type=Path)
    parser.add_argument('--kinetics-artifact-root', type=Path)
    args = parser.parse_args()
    export_visuals(args.output, args.gene_artifact_root, args.kinetics_artifact_root)
