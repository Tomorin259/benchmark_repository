"""Retain the accepted paired PCA design; export at 600 dpi by default."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

METHODS = ('veloae', 'scvelo_stochastic', 'TIvelo')
NAMES = ('veloAE', 'scVelo_sto', 'TIvelo')
COLORS = ('#C43C39', '#E17C05', '#B279A2')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-csv', type=Path, default=Path(__file__).resolve().parent / 'inputs/pca_comparison.csv')
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--dpi', type=int, default=600)
    args = parser.parse_args()
    with args.input_csv.open(newline='') as handle:
        rows = list(csv.DictReader(handle))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    width_inches = 180 / 25.4
    fig = plt.figure(figsize=(width_inches, width_inches * 420 / 1060))
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set(xlim=(0, 1060), ylim=(420, 0))
    ax.axis('off')
    pixel_to_pt = width_inches * 72 / 1060

    def line(x, y, u, v, color='#333333', width=1, opacity=1):
        ax.plot((x, u), (y, v), color=color, lw=width*pixel_to_pt,
                alpha=opacity, solid_capstyle='butt')

    def text(x, y, label, anchor='center', size=16, rotation=0):
        ax.text(x, y, label, ha=anchor, va='center' if rotation else 'baseline', color='#262626',
                fontsize=size*pixel_to_pt, fontfamily='DejaVu Sans', rotation=rotation)

    def dot(x, y, color, hollow=False, radius=4):
        ax.add_patch(Circle((x, y), radius, facecolor='white' if hollow else color,
                            edgecolor=color if hollow else '#333333', lw=pixel_to_pt))

    records = []
    for j, matrix in enumerate(('Ms', 'velocity')):
        left, top, bottom = 70 + j*530, 40, 325
        right = left + 435
        y = lambda q: bottom - float(q)/1.03*(bottom-top)
        for tick in (0, .2, .4, .6, .8, 1):
            line(left, y(tick), right, y(tick), '#E7E7E7', .7)
            text(left-10, y(tick)+5, f'{tick:.1f}', 'right', 14)
        line(left, top, left, bottom, width=1.4)
        line(left, bottom, right, bottom, width=1.4)
        text((left+right)/2, 23, 'Ms' if matrix == 'Ms' else 'Velocity', size=18)
        for i, (method, name, color) in enumerate(zip(METHODS, NAMES, COLORS)):
            subset = sorted((r for r in rows if r['method'] == method and
                             r['matrix_key'] == matrix), key=lambda r: r['dataset'])
            if len(subset) != 14 or len({r['dataset'] for r in subset}) != 14:
                raise ValueError(f'{matrix}/{method}: expected 14 unique datasets')
            center = left + 75 + i*145
            for k, row in enumerate(subset):
                raw, standardized = float(row['raw_pc1']), float(row['standardized_pc1'])
                if not (0 <= raw <= 1 and 0 <= standardized <= 1):
                    raise ValueError(f'Invalid explained variance: {row}')
                a, b = center-24+(k-6.5)*1.4, center+24+(k-6.5)*1.4
                line(a, y(raw), b, y(standardized), color, .8, .35)
                dot(a, y(raw), color, radius=3.7)
                dot(b, y(standardized), color, True, 3.7)
            for xx, key in ((center-24, 'raw_pc1'), (center+24, 'standardized_pc1')):
                mean = sum(float(r[key]) for r in subset)/14
                line(xx-12, y(mean), xx+12, y(mean), '#262626', 2)
            text(center, bottom+26, name, size=14)
            records.extend(subset)
    text(18, 185, 'PC1 explained variance', rotation=90)
    dot(390, 390, '#666666'); text(403, 395, 'Raw', 'left', 14)
    dot(485, 390, '#666666', True); text(498, 395, 'Standardized', 'left', 14)
    for suffix in ('png', 'pdf', 'svg'):
        fig.savefig(args.output_dir / f'pca_summary.{suffix}', dpi=args.dpi, facecolor='white')
    plt.close(fig)
    provenance = {'input_name': args.input_csv.name,
                  'input_sha256': hashlib.sha256(args.input_csv.read_bytes()).hexdigest(),
                  'dpi': args.dpi, 'width_mm': 180, 'paired_rows': len(records)}
    (args.output_dir / 'pca_summary_manifest.json').write_text(json.dumps(provenance, indent=2)+'\n')


if __name__ == '__main__':
    main()
