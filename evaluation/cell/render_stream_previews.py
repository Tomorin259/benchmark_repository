#!/usr/bin/env python3
"""Materialize static Cell previews from the registered, read-only outputs."""
from __future__ import annotations

import argparse
import gc
import hashlib
import json
from pathlib import Path

import h5py
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scvelo as scv

from prepare_cell_metrics import (align_cells_to_reference, align_genes_to_reference,
    apply_reference_graph, clear_derived_graphs,
    graphvelo_dynamo_parent_fields, inject_reference_state, load_reference_moments,
    prepare_velocity, project_forebrain_pca, read_inputs, resolve_source, sha256)
from reference_embeddings import array_sha256, load_forebrain_frames

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def materialize(artifact_root, moments_root, cache_root, output, datasets, manifest_root, methods=None):
    inputs = read_inputs(HERE / 'inputs.tsv')
    config = json.loads((HERE / 'config.json').read_text())
    palettes = json.loads((HERE / 'stream_preview_config.json').read_text())
    entries = []
    scv.settings.verbosity = 0
    for dataset, group in inputs.groupby('evaluation_dataset', sort=True):
        if datasets and dataset not in datasets:
            continue
        entries = []
        settings = palettes[dataset]
        reference_row = group[group.operational_id == 'scvelo_deterministic'].iloc[0]
        reference_path = resolve_source(artifact_root, reference_row.result_path)
        if sha256(reference_path) != reference_row.sha256:
            raise ValueError(f'Reference source changed: {dataset}')
        with h5py.File(reference_path, 'r') as reference:
            obs = reference['obs']
            ref_names = obs[obs.attrs['_index']].asstr()[:]
            labels = obs['clusters']
            label_values = labels['categories'].asstr()[:][labels['codes'][:]]
            frame = reference['obsm']['X_umap'][:, :2]
        if set(label_values) - set(settings['categories']):
            raise ValueError(f'Unexpected cell labels: {dataset}')
        if dataset == 'forebrain':
            frames = load_forebrain_frames()
            if not np.array_equal(ref_names, frames['obs_names']):
                raise ValueError('Forebrain reference and frozen PCA cell orders differ')
            frame = frames['X_pca'][:, :2]
        frame_contract = {
            'space': 'pca' if dataset == 'forebrain' else 'umap',
            'projection_dimensions': 30 if dataset == 'forebrain' else 2,
            'display_components': [1, 2],
            'coordinate_sha256': frames['metadata']['pca_coordinate_sha256']
                if dataset == 'forebrain' else array_sha256(frame),
        }
        protocol = {
            'renderer_sha256': sha256(Path(__file__)),
            'preparation_sha256': sha256(HERE / 'prepare_cell_metrics.py'),
            'frame_loader_sha256': sha256(HERE / 'reference_embeddings.py'),
            'config_sha256': sha256(HERE / 'config.json'),
            'palette_sha256': sha256(HERE / 'stream_preview_config.json'),
            'reference_source_sha256': reference_row.sha256,
            'scvelo_version': scv.__version__, 'frame': frame_contract,
        }
        if dataset == 'forebrain':
            protocol['reference_neighbors_sha256'] = sha256(HERE / 'reference_neighbors/forebrain.npz')
            protocol['reference_moments_sha256'] = sha256(moments_root / 'forebrain.npz')
            protocol['direct_pca_projection'] = False
        protocol_sha = hashlib.sha256(json.dumps(protocol, sort_keys=True).encode()).hexdigest()
        manifest = {'schema_version': 1, 'inputs_sha256': sha256(HERE / 'inputs.tsv'),
            'renderer_sha256': sha256(Path(__file__)),
            'palette_sha256': sha256(HERE / 'stream_preview_config.json'),
            'scvelo_version': scv.__version__, 'dataset_frames': {reference_row.dataset_id: frame_contract},
            'render_protocols': {reference_row.dataset_id: protocol}, 'entries': entries}
        fragment = manifest_root / f'{dataset}.json'
        previous = json.loads(fragment.read_text())['entries'] if fragment.exists() else []
        cached = {(e['dataset'],e['method']): e for e in previous}
        legend = [{'label': label, 'color': color} for label, color in
                  zip(settings['categories'], settings['palette'])]
        for _, row in group.iterrows():
            if methods and row.operational_id not in methods:
                continue
            source = resolve_source(artifact_root, row.result_path)
            prior = cached.get((row.dataset_id,row.config_id))
            if prior and prior.get('render_protocol_sha256') == protocol_sha and prior.get('frame') == frame_contract and prior['result_sha256'] == row.sha256 and source.stat().st_size == int(row.file_size_bytes) and source.stat().st_mtime <= fragment.stat().st_mtime and sha256(ROOT / prior['path']) == prior['sha256']:
                entries.append(prior)
                print(f'{len(entries)}/20 {row.dataset_id}/{row.config_id} resumed verified preview',flush=True)
                continue
            source_hash = sha256(source)
            if source_hash != row.sha256:
                raise ValueError(f'Method source changed: {dataset}/{row.operational_id}')
            existing = HERE / 'streams' / dataset / f'{row.operational_id}.png'
            provenance = 'retained_paper_panel'
            if existing.is_file() and row.operational_id != 'multivelo':
                image = existing
            else:
                image = output / row.dataset_id / f'{row.config_id}.png'
                image.parent.mkdir(parents=True, exist_ok=True)
                with h5py.File(source, 'r') as handle:
                    obs = handle['obs']
                    method_names = obs[obs.attrs['_index']].asstr()[:]
                positions = pd.Index(method_names).get_indexer(ref_names)
                reference_positions = np.flatnonzero(positions >= 0)
                method_positions = positions[reference_positions]
                if len(reference_positions) != int(row.n_common):
                    raise ValueError(f'Cell coverage differs: {dataset}/{row.operational_id}')
                xy = frame[reference_positions]
                if dataset != 'forebrain':
                    with h5py.File(source,'r') as handle:
                        native_xy = handle['obsm'][row.source_embedding_key][:, :2][method_positions]
                        velocity = handle['obsm'][row.velocity_embedding_key][:, :2][method_positions]
                    if not np.allclose(native_xy, xy, rtol=0, atol=1e-5):
                        raise ValueError(f'Noncanonical display frame: {dataset}/{row.operational_id}')
                    provenance = 'saved_native_projection_in_reference_frame'
                else:
                    method = ad.read_h5ad(source)
                    prepare_velocity(method, row.velocity_layer_key)
                    moments = load_reference_moments(dataset, moments_root)
                    method, cell_positions, _ = align_cells_to_reference(method, moments['obs_names'], row)
                    if row.operational_id != 'veloae':
                        extra = graphvelo_dynamo_parent_fields(method, dataset, inputs, artifact_root) \
                            if row.operational_id == 'graphvelo_dynamo' else None
                        method, gene_positions, _ = align_genes_to_reference(method, moments['var_names'], extra)
                        inject_reference_state(method, moments, cell_positions, gene_positions)
                    apply_reference_graph(method, dataset, cell_positions, HERE / 'reference_neighbors')
                    clear_derived_graphs(method)
                    scv.tl.velocity_graph(method, vkey='velocity', xkey='Ms', n_neighbors=30,
                        n_recurse_neighbors=1, mode_neighbors='distances', n_jobs=8)
                    project_forebrain_pca(method)
                    xy = np.asarray(method.obsm['X_pca'])[:, :2]
                    velocity = np.asarray(method.obsm['velocity_pca'])[:, :2]
                    provenance = 'fixed_reference30_pca30_projection_pc1_pc2'
                    del moments
                velocity = np.nan_to_num(velocity, nan=0, posinf=0, neginf=0)
                plot = ad.AnnData(X=np.zeros((len(xy), 1), dtype=np.float32))
                plot.obs['clusters'] = pd.Categorical(label_values[reference_positions],
                    categories=settings['categories'], ordered=True)
                plot.uns['clusters_colors'] = settings['palette']
                plot.obsm['X_umap'] = xy
                fig, ax = plt.subplots(figsize=(6, 4.5))
                scv.pl.velocity_embedding_stream(plot, basis='umap', X=xy, V=velocity,
                    color='clusters', palette=settings['palette'], title='', legend_loc='none',
                    frameon=False, show=False, ax=ax)
                # The reference bounds are fixed across all methods for a dataset.
                spread = np.ptp(frame, axis=0)
                ax.set_xlim(frame[:, 0].min() - spread[0]*.05, frame[:, 0].max() + spread[0]*.05)
                ax.set_ylim(frame[:, 1].min() - spread[1]*.05, frame[:, 1].max() + spread[1]*.05)
                fig.savefig(image, dpi=180, facecolor='white')
                plt.close(fig)
                if dataset == 'forebrain':
                    del method
                del plot
                gc.collect()
            entries.append({'dataset': row.dataset_id, 'method': row.config_id,
                'path': str(image.relative_to(ROOT)), 'sha256': sha256(image),
                'result_path': row.result_path, 'result_sha256': source_hash,
                'n_reference': int(row.n_reference), 'n_common': int(row.n_common),
                'projection': provenance, 'legend': legend, 'frame': frame_contract,
                'render_protocol_sha256': protocol_sha})
            (manifest_root / f'{dataset}.json').write_text(json.dumps(manifest, indent=2))
            print(f'{len(entries)}/280 {row.dataset_id}/{row.config_id} {provenance}', flush=True)

        manifest['renderer_sha256'] = sha256(Path(__file__))
        (manifest_root / f'{dataset}.json').write_text(json.dumps(manifest,indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--artifact-root', type=Path, required=True)
    parser.add_argument('--reference-moments-root', type=Path, required=True)
    parser.add_argument('--cache-root', type=Path, required=True)
    parser.add_argument('--datasets', nargs='+')
    parser.add_argument('--method', action='append', help='Select operational method IDs')
    parser.add_argument('--manifest-root', type=Path, required=True)
    args = parser.parse_args()
    args.manifest_root.mkdir(parents=True,exist_ok=True)
    materialize(args.artifact_root, args.reference_moments_root, args.cache_root,
                HERE / 'streams' / 'previews', args.datasets, args.manifest_root, args.method)

    if args.datasets is None and args.method is None:
        fragments = [json.loads(f.read_text()) for f in sorted(args.manifest_root.glob('*.json'))]
        combined = {k:v for k,v in fragments[0].items() if k!='entries'}
        combined['entries'] = [entry for fragment in fragments for entry in fragment['entries']]
        for key in ('dataset_frames', 'render_protocols'):
            combined[key] = {name: value for fragment in fragments
                for name, value in fragment[key].items()}
        if len(combined['entries']) != 280:
            raise ValueError('Expected all 280 registered Cell previews')
        (HERE / 'streams/preview_manifest.json').write_text(json.dumps(combined,indent=2))
