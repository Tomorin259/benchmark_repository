"""Align the shared noise-free ODE velocity reference by cell and gene name."""
from pathlib import Path
import csv
import hashlib
import json

import numpy as np


DEFAULT_TRUTH_ROOT = Path(__file__).resolve().parent / "truth"


def load_velocity_reference(scenario, cell_ids, gene_names, truth_root=None):
    path = Path(truth_root or DEFAULT_TRUTH_ROOT) / f"{scenario}_truth.npz"
    requested_cells = np.asarray(cell_ids).astype(str)
    requested_genes = np.asarray(gene_names).astype(str)
    with np.load(path, allow_pickle=False) as saved:
        cells = saved["cell_ids"].astype(str)
        genes = saved["gene_names"].astype(str)
        s = saved["velocity_s_true"].astype(np.float64)
        u = saved["velocity_u_true"].astype(np.float64)
    for label, names in [("reference cells", cells), ("reference genes", genes),
                         ("requested cells", requested_cells), ("requested genes", requested_genes)]:
        if names.ndim != 1 or len(set(names)) != len(names):
            raise ValueError(f"Nonunique or invalid {label}")
    if s.shape != (len(cells), len(genes)) or u.shape != s.shape:
        raise ValueError(f"Reference shape differs from identifiers: {path}")
    if not np.isfinite(s).all() or not np.isfinite(u).all():
        raise ValueError(f"Nonfinite reference velocities: {path}")
    ci, gi = {name: i for i, name in enumerate(cells)}, {name: i for i, name in enumerate(genes)}
    if not set(requested_cells).issubset(ci) or not set(requested_genes).issubset(gi):
        raise ValueError(f"Reference does not cover requested cells/genes: {path}; use the matching reference scope")
    indices = np.ix_([ci[name] for name in requested_cells], [gi[name] for name in requested_genes])
    return s[indices], u[indices]


def digest(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def materialize_clean_reference(repository_root, artifact_root, output_root, all_genes=False):
    """Export named ODE reference arrays from registered clean-state H5AD inputs."""
    import anndata as ad
    from evaluation.gene.config import SIM_CONFIGS

    repository_root, artifact_root, output_root = map(Path, (repository_root, artifact_root, output_root))
    with (repository_root / 'data/input_file_manifest.tsv').open() as handle:
        registered = {row['path']: row['sha256'] for row in csv.DictReader(handle, delimiter='\t')}
    if output_root.exists() and any(output_root.iterdir()):
        raise ValueError('Reference output directory must be empty')
    output_root.mkdir(parents=True, exist_ok=True)
    records = []
    for scenario in ['murk', 'branching', 'all']:
        relative = SIM_CONFIGS[scenario]['reference_source']
        source = artifact_root / relative
        source_sha = digest(source)
        if source_sha != registered.get(relative):
            raise ValueError('Registered clean-state input hash differs for ' + scenario)
        data = ad.read_h5ad(source, backed='r')
        try:
            required = {'clean_spliced', 'clean_unspliced', 'velocity_s_true', 'velocity_u_true',
                        'alpha_true', 'beta_true', 'gamma_true'}
            if not required.issubset(data.layers):
                raise ValueError('Input lacks clean-state ODE layers: ' + scenario)
            cells = np.asarray(data.obs_names).astype(str)
            genes = np.asarray(data.var_names if all_genes else SIM_CONFIGS[scenario]['genes']).astype(str)
            if len(set(cells)) != len(cells) or len(set(genes)) != len(genes):
                raise ValueError('Duplicate cell or gene identifier in ' + scenario)
            columns = data.var_names.get_indexer(genes)
            if np.any(columns < 0):
                raise ValueError('Input does not cover selected genes: ' + scenario)
            arrays = {key: np.asarray(data.layers[layer][:, columns], dtype=np.float64)
                      for key, layer in [('clean_s', 'clean_spliced'), ('clean_u', 'clean_unspliced'),
                                         ('velocity_s_true', 'velocity_s_true'), ('velocity_u_true', 'velocity_u_true'),
                                         ('alpha', 'alpha_true'), ('beta', 'beta_true'), ('gamma', 'gamma_true')]}
            if any(value.shape != (len(cells), len(genes)) or not np.isfinite(value).all()
                   for value in arrays.values()):
                raise ValueError('Clean-state array shape or finite-value contract differs: ' + scenario)
            if not (np.array_equal(arrays['velocity_s_true'], arrays['beta'] * arrays['clean_u'] - arrays['gamma'] * arrays['clean_s'])
                    and np.array_equal(arrays['velocity_u_true'], arrays['alpha'] - arrays['beta'] * arrays['clean_u'])):
                raise ValueError('Clean-state ODE layers differ for ' + scenario)
            output = output_root / (scenario + '_truth.npz')
            np.savez_compressed(output, cell_ids=cells, gene_names=genes,
                                **{k: arrays[k] for k in ['clean_u', 'clean_s', 'velocity_u_true', 'velocity_s_true']})
            records.append(dict(scenario=scenario, source=relative, source_sha256=source_sha,
                                output=output.name, output_sha256=digest(output),
                                shape=[len(cells), len(genes)], input_contract='clean_state_layers'))
        finally:
            data.file.close()
    manifest = dict(schema_version=1, truth='noise-free cellDancer ODE state',
                    generator_sha256=digest(__file__), coordinate_order=['spliced', 'unspliced'],
                    scenario_records=records)
    (output_root / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    return records


if __name__ == '__main__':
    import argparse
    import sys
    root = Path(__file__).resolve().parents[3]
    sys.path.insert(0, str(root))
    parser = argparse.ArgumentParser(description='Export registered clean-state simulation reference arrays.')
    parser.add_argument('--artifact-root', type=Path, required=True)
    parser.add_argument('--output-root', type=Path, required=True)
    parser.add_argument('--all-genes', action='store_true')
    args = parser.parse_args()
    materialize_clean_reference(root, args.artifact_root, args.output_root, args.all_genes)
