import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path

SCENARIOS = [('murk', 'Normal + MURK'), ('branching', 'Normal + branching'), ('all', 'Hybrid')]
METRICS = [('consistency', 'Consistency Score'), ('cbdir', 'CBDir'), ('iccoh', 'ICCoh'),
           ('transcosine', 'TransCosine'), ('batc', 'BATC')]


def main():
    root = Path(__file__).resolve().parent
    repository = root.parents[2]
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-root', type=Path, default=root / 'results')
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    output = args.output_root
    source = root / 'results/cell_metrics.csv'
    registry = repository / 'methods/configurations.tsv'
    with source.open(newline='') as handle:
        source_rows = list(csv.DictReader(handle))
    values = {(row['dataset'], row['method']): row for row in source_rows}
    with registry.open() as handle:
        names = {row['operational_id']: row['manuscript_name'] for row in csv.DictReader(handle, delimiter='\t')}
    with (root / 'inputs.tsv').open(newline='') as handle:
        input_methods = {row['operational_id'] for row in csv.DictReader(handle, delimiter='\t')}
    methods = [method for method in names if method in input_methods]
    expected_pairs = {(dataset, method) for dataset, _ in SCENARIOS for method in methods}
    if len(methods) != 9 or len(source_rows) != 27 or set(values) != expected_pairs:
        raise ValueError('Expected three scenarios, nine methods and five metrics')
    formatted = []
    blocks = []
    for method in methods:
        block = []
        for dataset, scenario in SCENARIOS:
            row = values[(dataset, method)]
            numeric = [float(row[key]) for key, _ in METRICS]
            if not all(math.isfinite(value) and -1 <= value <= 1 for value in numeric):
                raise ValueError('Expected finite cell-level scores within [-1, 1]')
            scores = [format(value, '.4f') for value in numeric]
            formatted.append([names[method], scenario] + scores)
            block.append([names[method], scenario] + scores)
        blocks.append(block)
    headers = ['Method', 'Dataset'] + [name for _, name in METRICS]
    csv_buffer = io.StringIO(newline='')
    writer = csv.writer(csv_buffer, lineterminator='\n')
    writer.writerow(headers)
    writer.writerows(formatted)
    lines = [r'\begin{table}[ht]', r'    \centering', r'    \resizebox{.75\textwidth}{!}{%',
             r'\small', r'\setlength{\tabcolsep}{4pt}', r'\begin{tabular}{llrrrrr}', r'\toprule',
             r'\multirow{2}{*}{Method} & \multirow{2}{*}{Dataset} & \multicolumn{5}{c}{Cell-level metrics} \\',
             r'\cmidrule(lr){3-7}',
             '& & ' + ' & '.join(name for _, name in METRICS) + r' \\', r'\midrule']
    for index, block in enumerate(blocks):
        if index:
            lines.append(r'\midrule')
        for row in block:
            lines.append(' & '.join([row[0].replace('_', r'\_')] + row[1:]) + r' \\')
    lines += [r'\bottomrule', r'\end{tabular}', r'    }',
              r'    \caption{\textbf{\normalsize Cell-level benchmarking results for nine RNA velocity method configurations across three simulated datasets.}\label{tab:suppSimulationCellMetrics}',
              r'    \footnotesize', r'    \normalfont',
              '    For each method and simulated dataset, we report Consistency Score, cross-boundary direction correctness (CBDir), in-cluster coherence (ICCoh), transition cosine similarity (TransCosine), and branching-aware trajectory consistency (BATC). The datasets are Normal + MURK, Normal + branching, and Hybrid. Higher values indicate better performance for all five metrics. Scores are reported to four decimal places.',
              r'    }', r'\end{table}']
    rounded_two = sum(format(float(values[(dataset, method)][key]), '.2f') == '1.00'
                      for dataset, _ in SCENARIOS for method in methods for key, _ in METRICS)
    rounded_four = sum(row[2:].count('1.0000') for row in formatted)
    manifest = {'source': str(source.relative_to(repository)), 'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
                'scenarios': [name for _, name in SCENARIOS], 'methods': methods, 'score_count': 135,
                'display_format': '.4f', 'source_precision': 'Unchanged full-precision CSV',
                'scores_rounding_to_one_at_two_decimals': rounded_two,
                'scores_rounding_to_one_at_four_decimals': rounded_four,
                'latex_dependencies': ['booktabs', 'multirow', 'graphicx'],
                'table_label': 'tab:suppSimulationCellMetrics',
                'style_source': 'Online merged.tex, Supplemental Tables, tab:tableThreeA--E',
                'column_order': headers, 'row_grouping': 'Method; three datasets per method',
                'encoding': 'Raw scores; inline LaTeX table; no ranks, heatmap, table image or PDF export'}
    exports = {'simulation_cell_scores_4dp.csv': csv_buffer.getvalue(),
               'simulation_cell_scores.tex': '\n'.join(lines) + '\n',
               'table_manifest.json': json.dumps(manifest, indent=2) + '\n'}
    if args.check:
        for filename, content in exports.items():
            path = output / filename
            if not path.is_file() or path.read_bytes() != content.encode('utf-8'):
                raise ValueError(f'Stale supplementary table export: {path}')
    else:
        output.mkdir(parents=True, exist_ok=True)
        for filename, content in exports.items():
            (output / filename).write_text(content, encoding='utf-8')
    print('Table: 135 raw scores displayed to four decimals; {} -> {} scores round to one'.format(rounded_two, rounded_four))


if __name__ == '__main__':
    main()
