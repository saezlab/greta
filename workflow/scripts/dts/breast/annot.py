import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import sys
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-a','--path_tmp', required=True)
parser.add_argument('-b','--path_gid', required=True)
parser.add_argument('-c','--path_rna', required=True)
parser.add_argument('-d','--path_ann', required=True)
parser.add_argument('-e','--sample_ids', required=True, nargs='+')
args = vars(parser.parse_args())

path_tmp = args['path_tmp']
path_gid = args['path_gid']
path_rna = args['path_rna']
path_ann = args['path_ann']
sample_ids = args['sample_ids']

# Read
rna = ad.read_h5ad(path_tmp)
geneids = pd.read_csv(path_gid).set_index('id')['symbol'].to_dict()
# Format obs
rna.obs['ID'] = 'ID' + rna.obs.index.str.split('_').str[1]
rna.obs = rna.obs[['ID', 'celltype']].rename(columns={'ID': 'batch'})
rna.obs.index = [b + '_' + i.split('-')[0].split('_')[-1] for i, b in zip(rna.obs.index, rna.obs['batch'])]
# Basic filter
sc.pp.filter_cells(rna, min_genes=100)
sc.pp.filter_genes(rna, min_cells=3)
# Filter faulty gene symbols
ensmbls = np.array([geneids[g] if g in geneids else '' for g in rna.var_names])
msk = ensmbls != ''
rna = rna[:, msk].copy()
rna.var['symbol'] = ensmbls[msk]
# Remove duplicated genes based on num of cells
to_remove = []
for dup in rna.var['symbol'].values[rna.var['symbol'].duplicated()]:
    tmp = rna.var[rna.var['symbol'] == dup]
    max_idx = tmp['n_cells'].idxmax()
    to_remove.extend(tmp.index[tmp.index != max_idx].values)
rna = rna[:, ~rna.var_names.isin(to_remove)].copy()
rna.var_names = rna.var['symbol']
# Clean
del rna.obs['n_genes']
del rna.var
del rna.uns
del rna.obsm
# Filter low abundant cells and sample ids
msk_ids = rna.obs['batch'].isin(sample_ids)
tobs = rna.obs[msk_ids]
ncells = tobs.groupby('celltype').size()
ctypes = set(ncells.index[ncells >= 50])
msk_obs = (rna.obs['celltype'].isin(ctypes)) & msk_ids
rna = rna[msk_obs, :].copy()

# Write
rna.obs.to_csv(path_ann)
del rna.obs
rna.write(path_rna)
