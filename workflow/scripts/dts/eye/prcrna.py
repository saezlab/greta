import scanpy as sc
import numpy as np
import pandas as pd
import os
import sys

path_geneids = sys.argv[1]
path_ann = sys.argv[2]
path_out = sys.argv[3]

# Read
geneids = pd.read_csv(path_geneids).set_index('symbol')['id'].to_dict()
path_dir = os.path.dirname(path_out)
rna = sc.read_10x_mtx(path_dir, prefix='GSM5866081_')
obs = pd.read_csv(path_ann, index_col=0)
bar2sam = dict()
for i in obs.index:
    sam, bar = i.split('_')
    bar2sam[bar] = sam
# Filter by obs
rna.obs_names = [b.split('_')[1].split('-')[0] for b in rna.obs_names]
msk = rna.obs_names.isin(bar2sam)
rna = rna[msk, :].copy()
rna.obs_names = [f'{bar2sam[b]}_{b}' for b in rna.obs_names]
rna = rna[obs.index, :].copy()
# Filter faulty gene symbols
ensmbls = np.array([geneids[g] if g in geneids else '' for g in rna.var_names])
msk = ensmbls != ''
rna = rna[:, msk].copy()
rna.var['gene_ids'] = ensmbls[msk]
# Basic QC
sc.pp.filter_cells(rna, min_genes=100)
sc.pp.filter_genes(rna, min_cells=3)
del rna.obs['n_genes']
# Remove duplicated genes based on num of cells
to_remove = []
for dup in rna.var.index[rna.var.index.duplicated()]:
    tmp = rna.var.loc[dup]
    max_idx = tmp.set_index('gene_ids')['n_cells'].idxmax()
    to_remove.extend(tmp['gene_ids'][tmp['gene_ids'] != max_idx].values)
rna = rna[:, ~rna.var['gene_ids'].isin(to_remove)].copy()
del rna.obs
del rna.var
# Write
rna.write(path_out)
