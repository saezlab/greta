import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--path_atac', required=True)
parser.add_argument('-o','--path_annot', required=True)
args = vars(parser.parse_args())


path_atac = args['path_atac']
path_annot = args['path_annot']

# Get obs
atac = ad.read_h5ad(path_atac).obs
atac = atac[atac['region'] == 'LV'][['combinedID', 'cell_type']].copy()
atac[['sangerID', 'batch']] = atac['combinedID'].str.split('_', expand=True)
atac.index = [b + '_' + i.split('-')[0].split('_')[-1] for i, b in zip(atac.index, atac['batch'])]
atac = atac.rename(columns={'cell_type': 'celltype', 'sangerID': 'sangerid'})
atac = atac[['celltype', 'batch', 'sangerid']]
ctype_counts = atac.groupby('celltype', as_index=False).size()
ctypes = ctype_counts[ctype_counts['size'] >= 100]['celltype'].values.astype(str)
atac = atac[atac['celltype'].isin(ctypes)]

# Write
atac.to_csv(path_annot)