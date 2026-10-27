#!/usr/bin/env python3
"""Render the canonical five-metric Cell figures from split result tables."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

PLOT_RANDOM_SEED = 123
MODULE_DIR = Path(__file__).resolve().parent
DEFAULT_OUTDIR = MODULE_DIR / 'figures'
DEFAULT_METRICS_CSV = MODULE_DIR / 'results/cell_metrics.csv'
DEFAULT_SUPPLEMENTARY_METRICS_CSV = (
    MODULE_DIR / 'results/cell_metrics_supplementary.csv'
)
OUTDIR = DEFAULT_OUTDIR
METRICS_CSV = DEFAULT_METRICS_CSV
SUPPLEMENTARY_METRICS_CSV = DEFAULT_SUPPLEMENTARY_METRICS_CSV

METHOD_NAME_MAPPING = {
    'scvelo_deterministic': 'scVelo_det',
    'scvelo_dynamic': 'scVelo_dyn',
    'scvelo_stochastic': 'scVelo_sto',
    'unitvelo': 'UniTVelo',
    'unitvelo_mode2': 'UniTVelo(mode2)',
    'celldancer': 'cellDancer',
    'TFvelo': 'TFvelo',
    'dynamo': 'Dynamo(RNA_seq)',
    'deepvelo': 'DeepVelo',
    'SDEvelo': 'SDEvelo',
    'DeepKINET': 'DeepKINET',
    'cell2fate': 'cell2fate',
    'TIvelo': 'TIvelo',
    'veloVI': 'veloVI',
    'velovi': 'veloVI',
    'multivelo': 'MultiVelo (RNA_only)',
    'deepvelo2022': 'DeepVelo_scv',
    'graphvelo_scv': 'GraphVelo_scv',
    'graphvelo_dynamo': 'GraphVelo_dyn',
    'TSvelo': 'TSvelo',
    'veloae': 'veloAE',
}

METHOD_ORDER = [
    'scVelo_det', 'scVelo_sto', 'scVelo_dyn', 'MultiVelo (RNA_only)',
    'UniTVelo', 'UniTVelo(mode2)', 'Dynamo(RNA_seq)', 'TIvelo', 'cell2fate',
    'TFvelo', 'veloAE', 'veloVI', 'DeepKINET', 'TSvelo', 'cellDancer',
    'DeepVelo', 'SDEvelo', 'DeepVelo_scv', 'GraphVelo_scv',
    'GraphVelo_dyn',
]

DATASET_NAME_MAPPING = {
    'gastrulation': 'Gastrulation',
    'DentateGyrus_lamanno': 'Dentategyrus_lamanno',
    'dentategyrus_lamanno': 'Dentategyrus_lamanno',
    'dentategyrus': 'Dentategyrus',
    'bonemarrow': 'Human bone marrow',
    'capo_pre': 'capofibro (VeloCycle)',
    'capofibro': 'capofibro (VeloCycle)',
    'retina': 'Mouse retina',
    'mousebrain': 'mousebrain',
    'forebrain': 'ForebrainGlut',
    'hematopoiesis': 'Human hematopoiesis',
    'pancreas': 'Pancreas',
    'neuron': 'scNT_neuron',
    'organoid': 'organoid',
    'rpe1': 'rpe1',
    'rpe1_gen': 'rpe1 (VeloCycle)',
}

BENCHMARK_COLORS = {
    'scVelo_det': '#1F77B4',
    'scVelo_sto': '#FF7F0E',
    'Dynamo(RNA_seq)': '#2CA02C',
    'veloAE': '#D62728',
    'scVelo_dyn': '#9467BD',
    'MultiVelo (RNA_only)': '#8C564B',
    'UniTVelo': '#E377C2',
    'UniTVelo(mode2)': '#7F7F7F',
    'veloVI': '#BCBD22',
    'cell2fate': '#17BECF',
    'SDEvelo': '#0B559F',
    'DeepKINET': '#FA8F12',
    'cellDancer': '#5BA658',
    'DeepVelo': '#C9383B',
    'TFvelo': '#806BA8',
    'TSvelo': '#A05D56',
    'TIvelo': '#D16FA3',
    'DeepVelo_scv': '#5F6B6D',
    'GraphVelo_scv': '#B5BD61',
    'GraphVelo_dyn': '#2E91E5',
}

METRICS = [
    {
        'key': 'consistency',
        'column': 'consistency',
        'title': 'Consistency',
        'source': 'primary',
        'xlim': (0.1, 1.05),
        'baseline': None,
    },
    {
        'key': 'iccoh',
        'column': 'iccoh',
        'title': 'ICCoh',
        'source': 'primary',
        'xlim': (0.2, 1.05),
        'baseline': None,
    },
    {
        'key': 'cbdir',
        'column': 'cbdir',
        'title': 'CBDir',
        'source': 'primary',
        'xlim': (-1.05, 1.05),
        'baseline': 0,
    },
    {
        'key': 'transcosine',
        'column': 'transcosine',
        'title': 'TransCosine',
        'source': 'primary',
        'xlim': (-0.03, 0.72),
        'baseline': None,
    },
    {
        'key': 'batc',
        'column': 'batc',
        'title': 'BATC',
        'source': 'primary',
        'xlim': (-1.05, 1.05),
        'baseline': 0,
    },
]

HIGHDIM_CBDIR = {
    'key': 'highdim_cbdir',
    'column': 'highdim_cbdir',
    'title': 'CBDir-HD',
    'source': 'supplementary',
    'xlim': (-0.7, 0.7),
    'baseline': 0,
}


def load_metric(metric):
    if metric['source'] == 'primary':
        data = pd.read_csv(METRICS_CSV)
    elif metric['source'] == 'supplementary':
        data = pd.read_csv(SUPPLEMENTARY_METRICS_CSV)
    else:
        raise ValueError(metric['source'])

    required = {'dataset', 'method', metric['column']}
    missing = required - set(data.columns)
    if missing:
        raise ValueError(f"Missing columns in {metric['source']}: {sorted(missing)}")

    df = data[['dataset', 'method', metric['column']]].rename(columns={metric['column']: 'score'}).copy()
    df = df.dropna(subset=['score'])
    df['method'] = df['method'].replace(METHOD_NAME_MAPPING)
    df['dataset'] = df['dataset'].replace(DATASET_NAME_MAPPING)
    order = df.groupby('method')['score'].mean().sort_values(ascending=False).index.tolist()
    means = df.groupby('method')['score'].mean()
    return df, order, means


def draw_metric(ax, metric, compact=False):
    df, order, means = load_metric(metric)
    sns.boxplot(
        y='method',
        x='score',
        data=df,
        palette=BENCHMARK_COLORS,
        order=order,
        width=0.42,
        linewidth=0.85 if compact else 1.05,
        showfliers=False,
        boxprops={'alpha': 0.62},
        whiskerprops={'linewidth': 0.8 if compact else 1.0},
        capprops={'linewidth': 0.8 if compact else 1.0},
        medianprops={'linewidth': 0.9 if compact else 1.1, 'color': 'black'},
        ax=ax,
    )

    np.random.seed(PLOT_RANDOM_SEED)
    sns.stripplot(
        y='method',
        x='score',
        data=df,
        hue='method',
        palette=BENCHMARK_COLORS,
        order=order,
        size=2.2 if compact else 2.8,
        jitter=0.18,
        alpha=0.86,
        linewidth=0.35 if compact else 0.45,
        edgecolor='black',
        legend=False,
        ax=ax,
    )

    for index, method in enumerate(order):
        ax.scatter(
            means[method],
            index,
            color='#FFD700',
            s=10 if compact else 18,
            marker='D',
            edgecolor='black',
            linewidth=0.35 if compact else 0.55,
            zorder=10,
        )

    if metric['baseline'] is not None:
        ax.axvline(metric['baseline'], color='gray', linestyle='-', linewidth=1.0, zorder=0, alpha=0.55)

    ax.set_title(metric['title'], fontsize=8.5 if compact else 10, fontweight='bold', pad=4)
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_xlim(*metric['xlim'])
    ax.tick_params(axis='y', labelsize=6.4 if compact else 7.6, length=0, pad=1)
    ax.tick_params(axis='x', labelsize=5.8 if compact else 6.8, length=2, pad=1)
    ax.grid(True, axis='x', linestyle='--', alpha=0.38, linewidth=0.45)
    ax.grid(False, axis='y')
    for spine in ['top', 'right', 'left']:
        ax.spines[spine].set_visible(False)
    ax.spines['bottom'].set_linewidth(0.7)


def save_individual_panels():
    sns.set_theme(style='ticks', rc={'font.family': 'DejaVu Sans'})
    for metric in METRICS:
        fig, ax = plt.subplots(figsize=(3.2, 5.2))
        draw_metric(ax, metric, compact=False)
        fig.subplots_adjust(left=0.34, right=0.985, top=0.925, bottom=0.075)
        for ext in ['pdf', 'svg', 'png']:
            fig.savefig(OUTDIR / f"fig2_cell_metric_{metric['key']}.{ext}", dpi=600, transparent=False)
        plt.close(fig)


def save_five_panel():
    sns.set_theme(style='ticks', rc={'font.family': 'DejaVu Sans'})
    fig, axes = plt.subplots(1, len(METRICS), figsize=(15.8, 5.2), sharey=False)
    for ax, metric in zip(axes, METRICS):
        draw_metric(ax, metric, compact=True)
    fig.subplots_adjust(left=0.075, right=0.995, top=0.925, bottom=0.075, wspace=0.62)
    for ext in ['pdf', 'svg', 'png']:
        fig.savefig(OUTDIR / f"fig2_cell_metrics_5panel.{ext}", dpi=600, transparent=False)
    plt.close(fig)


def save_highdim_cbdir():
    sns.set_theme(style='ticks', rc={'font.family': 'DejaVu Sans'})
    fig, ax = plt.subplots(figsize=(3.2, 5.2))
    draw_metric(ax, HIGHDIM_CBDIR, compact=False)
    fig.subplots_adjust(left=0.34, right=0.985, top=0.925, bottom=0.075)
    for ext in ['pdf', 'svg', 'png']:
        fig.savefig(OUTDIR / f"fig2_cell_metric_{HIGHDIM_CBDIR['key']}.{ext}", dpi=600, transparent=False)
    plt.close(fig)


def save_matched_panels():
    """Render S9 Cell panels from the seven matched-output comparisons."""
    global METRICS_CSV
    primary = pd.read_csv(METRICS_CSV)
    matched = pd.read_csv(MODULE_DIR / 'gene_matched/results/cell_metrics.csv')
    methods = {'scvelo_dynamic', 'deepvelo2022', 'graphvelo_scv', 'TIvelo',
               'celldancer', 'deepvelo', 'unitvelo', 'veloVI', 'multivelo'}
    selected = primary[primary.dataset.isin(matched.dataset) & primary.method.isin(methods)]
    selected = selected.set_index(['dataset', 'method'])
    selected.update(matched.set_index(['dataset', 'method']))
    if len(selected) != 63:
        raise ValueError('S9 requires nine methods on seven datasets')
    selected.reset_index().to_csv(OUTDIR / 'cell_scores.csv', index=False)
    source = METRICS_CSV
    METRICS_CSV = OUTDIR / 'cell_scores.csv'
    try:
        for i, metric in enumerate(METRICS, 1):
            fig, ax = plt.subplots(figsize=(3.2, 2.7))
            draw_metric(ax, metric, compact=True)
            fig.subplots_adjust(left=.42, right=.98, bottom=.14, top=.89)
            for suffix in ['png', 'pdf', 'svg']:
                fig.savefig(OUTDIR / f"cell_panel_{i:02d}_{metric['key']}.{suffix}",
                            dpi=600, facecolor='white')
            plt.close(fig)
    finally:
        METRICS_CSV = source


def main():
    global OUTDIR, METRICS_CSV, SUPPLEMENTARY_METRICS_CSV
    parser = argparse.ArgumentParser(
        description='Render the formal five-panel and CBDir-HD cell-metric figures.'
    )
    parser.add_argument('--metrics', type=Path, default=DEFAULT_METRICS_CSV)
    parser.add_argument(
        '--supplementary-metrics',
        type=Path,
        default=DEFAULT_SUPPLEMENTARY_METRICS_CSV,
    )
    parser.add_argument('--output', type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument('--gene-matched', action='store_true', help='Render the five S9 Cell panels')
    args = parser.parse_args()
    METRICS_CSV = args.metrics.resolve()
    SUPPLEMENTARY_METRICS_CSV = args.supplementary_metrics.resolve()
    OUTDIR = args.output.resolve()
    OUTDIR.mkdir(parents=True, exist_ok=True)

    if args.gene_matched:
        save_matched_panels()
    else:
        save_individual_panels()
        save_five_panel()
        save_highdim_cbdir()
    print(f"Wrote outputs to {OUTDIR}")
    for path in sorted(OUTDIR.glob('fig2_cell_metric*')):
        print(path)


if __name__ == '__main__':
    main()
