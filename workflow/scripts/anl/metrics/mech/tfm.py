import scanpy as sc
import pandas as pd
import numpy as np
import mudata as mu
import os
import h5py
import sys


path_mdata = sys.argv[1]
path_tfs = sys.argv[2]
path_out = sys.argv[3]

# Read
tfs = pd.read_csv(path_tfs, header=None)[0].values
rna = mu.read(os.path.join(path_mdata, 'mod', 'rna'))

# Filter and update
inter = rna.var_names.intersection(tfs)
rna = rna[:, inter].copy()
rna.obs = mu.read(path_mdata).obs.loc[:, ['celltype']].copy()

# Extract DEG tfs
sc.tl.rank_genes_groups(rna, groupby='celltype', method='wilcoxon')
df = sc.get.rank_genes_groups_df(rna, group=None)

# Filter results
df = df[(df['pvals_adj'] < 2.22e-16) & (df['logfoldchanges'] > 2.)]
n_group = df.groupby('group', as_index=False).size()
n_group = n_group[n_group['size'] >= 1]
groups = n_group['group'].values
df['group'] = df['group'].astype(str)
df = df[df['group'].isin(groups)]
df = df[['group', 'names']]
df.columns = ['celltype', 'tf']

# Write
df.to_csv(path_out, index=False)
