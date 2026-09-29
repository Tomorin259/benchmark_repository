#!/usr/bin/env python3
"""Execute the pinned Barile notebook on an explicit Gastrulation input.

Only the input read is redirected; numerical notebook cells remain unchanged.
This is a reconstruction check, not a replacement of the frozen gene list.
"""
import argparse
import hashlib
import json
import os
import platform
from pathlib import Path

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--input', type=Path, required=True)
parser.add_argument('--input-id', required=True, help='Portable artifact locator')
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
source = HERE / 'upstream/search_murk_mouse_fig_3d_and_rank.ipynb'
input_file = args.input.resolve()
input_path = str(input_file)
output = args.output.resolve()
output.mkdir(parents=True, exist_ok=False)
notebook = json.loads(source.read_text())
os.chdir(str(output))
namespace = {'input_path': input_path}
for index, cell in enumerate(notebook['cells']):
    if cell['cell_type'] != 'code':
        continue
    code = ''.join(cell['source'])
    code = code.replace("dir_data + 'adata_umap_pca.h5'", 'input_path')
    exec(compile(code, 'Barile notebook cell {}'.format(index), 'exec'), namespace)

def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()

import numpy, pandas, scipy, scanpy, scvelo, anndata
metadata = {
    'status': 'reconstruction_check_not_manuscript_reselection',
    'input_id': args.input_id,
    'input_sha256': sha256(input_file),
    'upstream_commit': '1e0de6338a85a559478627927c24bdc72d2d8919',
    'notebook_sha256': sha256(source),
    'input_redirection_only': True,
    'versions': {m.__name__: m.__version__ for m in
                 [numpy, pandas, scipy, scanpy, scvelo, anndata]},
    'python': platform.python_version(),
    'pythonhashseed': os.environ.get('PYTHONHASHSEED'),
    'n_candidates': len(namespace['df_fin']),
}
(output / 'manifest.json').write_text(json.dumps(metadata, indent=2) + '\n')
