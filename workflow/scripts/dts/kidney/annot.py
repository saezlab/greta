import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
import anndata as ad
import sys
path_h5 = sys.argv[1]
path_obs = sys.argv[2]
path_adata = sys.argv[3]
# Read
with h5py.File(path_h5, 'r') as f:
    data = f['/assays/RNA/counts/data'][:]
    indices = f['/assays/RNA/counts/indices'][:]
    indptr = f['/assays/RNA/counts/indptr'][:]
    X = csr_matrix((data, indices, indptr))
    genes = f['assays']['RNA']['features'][:].astype('U')
    obs = pd.DataFrame(index=f['meta.data']['_index'][:].astype('U'))
    obs['batch'] = f['meta.data']['experiment'][:].astype('U')
    obs['celltype'] = f['meta.data']['subclass.l1'][:].astype('U')
# Format obs
obs['batch'] = [b.split('_')[1] for b in obs['batch']]
obs[['library', 'barcode']] = [i.split('_') for i in obs.index]
obs.index = [f'{batch}_{barcode}' for batch, barcode in zip(obs['batch'], obs['barcode'])]
obs.index = [i.replace('-1', '') for i in obs.index]
obs = obs.drop(columns=['library', 'barcode'])
# Filter low abundant cells
ncells = obs.groupby('celltype').size()
ctypes = set(ncells.index[ncells >= 100])
msk_obs = obs['celltype'].isin(ctypes)
X = X[msk_obs, :]
obs = obs.loc[msk_obs, :]
var = pd.DataFrame(index=genes)
adata = ad.AnnData(X=X, obs=obs.drop(columns=obs.columns), var=var)
# Write
obs.to_csv(path_obs)
adata.write(path_adata)
