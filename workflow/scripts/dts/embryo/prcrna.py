import h5py
import scipy.sparse as sps
import anndata as ad
import pandas as pd
import numpy as np
import os
import sys

path_h5 = sys.argv[1]
path_ann = sys.argv[2]
path_gid = sys.argv[3]

# Process annot
df = pd.read_csv(path_ann, index_col=0)
df = df[['cell_assignment']]
df['batch'] = [i.split('_')[0].upper() for i in df.index]
df.index = [b + '_' + i.split('_')[2].replace('-1', '') for b, i in zip(df['batch'], df.index)]
df = df.rename(columns={'cell_assignment': 'celltype'})
df = df[['batch', 'celltype']]
# Read anndata
with h5py.File(path_h5, "r") as f:
    obs = f['obs'][:].astype('U')
    var = f['var'][:].astype('U')
    X = sps.csc_matrix((f["x"][:], f["i"][:], f["p"][:])).tocsr().T
obs = pd.DataFrame(index=obs)
var = pd.DataFrame(index=var)
rna = ad.AnnData(X=X, obs=obs, var=var)
rna.obs_names = [i.replace('_double_', '_').upper().replace('-1', '') for i in rna.obs_names]
# Filter faulty gene symbols
geneids = pd.read_csv(path_gid).set_index('symbol')['id'].to_dict()
ensmbls = np.array([geneids[g] if g in geneids else '' for g in rna.var_names])
msk = ensmbls != ''
rna = rna[:, msk].copy()
rna.var['gene_ids'] = ensmbls[msk]
# Match ann
rna = rna[df.index, :].copy()
# Write
df.to_csv(path_ann)
path_rna = os.path.join(os.path.dirname(path_ann), 'rna.h5ad')
rna.write(path_rna)
