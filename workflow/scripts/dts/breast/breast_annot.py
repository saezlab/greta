import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--path_h5ad', required=True)
parser.add_argument('-o','--path_annot', required=True)
args = vars(parser.parse_args())

path_h5ad = args['path_h5ad']
path_annot = args['path_annot']

# Get obs
rna_obs = ad.read_h5ad(path_h5ad).obs

rna = rna_obs[['donor_id', 'Pool', 'cell_type']]
rna['ID'] = 'ID' + rna.index.str.split('_').str[1]
rna.index = [b + '_' + i.split('-')[0].split('_')[-1] for i, b in zip(rna.index, rna['ID'])]
rna = rna.rename(columns={'cell_type': 'celltype', 'ID': 'batch'})
rna = rna[['celltype', 'batch', 'donor_id', 'Pool']]
ctype_counts = rna.groupby('celltype', as_index=False).size()
ctypes = ctype_counts[ctype_counts['size'] >= 100]['celltype'].values.astype(str)
rna = rna[rna['celltype'].isin(ctypes)]

# Write
rna.to_csv(path_annot)