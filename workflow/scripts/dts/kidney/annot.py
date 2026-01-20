import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
import anndata as ad
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-a','--path_h5', required=True)
parser.add_argument('-b','--path_obs', required=True)
parser.add_argument('-c','--path_adata', required=True)
parser.add_argument('-d','--sample_ids', required=True, nargs='+')
args = vars(parser.parse_args())

path_h5 = args['path_h5']
path_obs = args['path_obs']
path_adata = args['path_adata']
sample_ids = args['sample_ids']

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
msk_ids = obs['batch'].isin(sample_ids)
tobs = obs[msk_ids]
ncells = tobs.groupby('celltype').size()
ctypes = set(ncells.index[ncells >= 50])
msk_obs = (obs['celltype'].isin(ctypes)) & msk_ids
X = X[msk_obs, :]
obs = obs.loc[msk_obs, :]
var = pd.DataFrame(index=genes)
adata = ad.AnnData(X=X, obs=obs.drop(columns=obs.columns), var=var)
# Write
obs.to_csv(path_obs)
adata.write(path_adata)
