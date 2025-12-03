#!/usr/bin/env python3
import argparse
import sys
import scanpy as sc
import pandas as pd

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-a', '--path_h5', required=True)
parser.add_argument('-b', '--path_annot', required=True)
parser.add_argument('-c', '--batch', required=True)
parser.add_argument('-d', '--path_out', required=True)
args = vars(parser.parse_args())

path_h5 = args['path_h5']
path_annot = args['path_annot']
batch = args['batch']
path_out = args['path_out']

adata = sc.read_10x_h5(path_h5)

if not adata.var_names.is_unique:
    adata.var_names_make_unique()

if 'feature_types' in adata.var:
    mask = adata.var['feature_types'] == 'Gene Expression'
    adata = adata[:, mask].copy()

new_barcodes = []
for bc in adata.obs_names:
    parts = bc.rsplit('-', 1)
    barcode_base = parts[0]
    new_barcodes.append(f'{batch}_{barcode_base}')

adata.obs_names = new_barcodes
print(f'[h5_to_h5ad] Converted {len(new_barcodes)} barcodes', file=sys.stderr)

annot = pd.read_csv(path_annot, index_col=0)
common = adata.obs_names.intersection(annot.index)

adata = adata[common, :].copy()
adata.obs['batch'] = annot.loc[adata.obs_names, 'batch']
adata.obs['celltype'] = annot.loc[adata.obs_names, 'celltype']

adata.write(path_out)
print(f'[h5_to_h5ad] Wrote {len(common)} cells to {path_out}', file=sys.stderr)
