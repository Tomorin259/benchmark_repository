#!/usr/bin/env python3
"""Compare full-precision numerical fields and saved row order, separately."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--historical', type=Path, required=True)
parser.add_argument('--rerun', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
old, new = pd.read_csv(args.historical), pd.read_csv(args.rerun)
old = old.rename(columns={old.columns[0]: 'gene'}).set_index('gene')
new = new.rename(columns={new.columns[0]: 'gene'}).set_index('gene')
if not old.index.is_unique or not new.index.is_unique:
    raise ValueError('Duplicate gene identifiers')
common = old.index.intersection(new.index)
columns = [c for c in old if c in new and pd.api.types.is_numeric_dtype(old[c])
           and pd.api.types.is_numeric_dtype(new[c])]
comparison = {}
for column in columns:
    a, b = old.loc[common, column].values, new.loc[common, column].values
    finite = np.isfinite(a) & np.isfinite(b)
    comparison[column] = {
        'n_common': len(common), 'n_both_finite': int(finite.sum()),
        'max_abs_difference': float(np.max(np.abs(a[finite] - b[finite]))) if finite.any() else None,
        'n_exact_equal': int(np.sum(a == b)),
        'n_abs_difference_le_1e_6': int(np.sum(np.isclose(a, b, rtol=0, atol=1e-6))),
        'n_abs_difference_le_1e_4': int(np.sum(np.isclose(a, b, rtol=0, atol=1e-4))),
    }
old_rank = pd.Series(np.arange(len(old)), index=old.index)
new_rank = pd.Series(np.arange(len(new)), index=new.index)
primary = 'scaled Δm' if 'scaled Δm' in columns else 'slope_gap' if 'slope_gap' in columns else 'separation_gap'
same_set = set(old.index) == set(new.index)
tie_equivalent = bool(same_set and primary in columns
                     and old[primary].is_monotonic_decreasing and new[primary].is_monotonic_decreasing
                     and np.array_equal(old.loc[common, primary].values, new.loc[common, primary].values))
result = {
    'historical_file': args.historical.name, 'rerun_file': args.rerun.name,
    'historical_sha256': hashlib.sha256(args.historical.read_bytes()).hexdigest(),
    'rerun_sha256': hashlib.sha256(args.rerun.read_bytes()).hexdigest(),
    'n_historical': len(old), 'n_rerun': len(new), 'n_common': len(common),
    'missing_from_rerun': sorted(set(old.index) - set(new.index)),
    'new_in_rerun': sorted(set(new.index) - set(old.index)),
    'row_order_exact': old.index.tolist() == new.index.tolist(),
    'primary_sort_column': primary,
    'row_order_differences_only_within_equal_score_ties': tie_equivalent,
    'tie_aware_primary_score_ranks_equal': bool(same_set and primary in columns and
        np.array_equal(old[primary].rank().loc[common].values, new[primary].rank().loc[common].values)),
    'common_gene_rank_spearman': float(spearmanr(old_rank.loc[common], new_rank.loc[common])[0]) if len(common) > 1 else None,
    'top50_overlap': len(set(old.index[:50]) & set(new.index[:50])),
    'numeric_columns': comparison,
}
args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
print(json.dumps({k: v for k, v in result.items() if k not in ['missing_from_rerun', 'new_in_rerun']}, indent=2))
