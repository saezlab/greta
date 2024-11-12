import numpy as np
import pandas as pd
import snapatac2 as snap
import scanpy as sc
import mudata as md
import scanpy.external as sce
from scipy.sparse import issparse
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--path_input', required=True)
parser.add_argument('-c','--celltypes', required=True)
parser.add_argument('-s','--n_sample', required=True)
parser.add_argument('-d','--seed', required=True)
parser.add_argument('-g','--n_hvg', required=True)
parser.add_argument('-r','--n_hvr', required=True)
parser.add_argument('-t','--root', required=True)
parser.add_argument('-o','--path_output', required=True)
args = vars(parser.parse_args())

path_input = args['path_input']
celltypes = args['celltypes']
n_sample = int(args['n_sample'])
seed = int(args['seed'])
n_hvg = int(args['n_hvg'])
n_hvr = int(args['n_hvr'])
root = args['root']
path_output = args['path_output']

# Read
mdata = md.read_h5mu(path_input)

# Filter celltypes
if celltypes != 'all':
    celltypes = celltypes.split(';')
    mdata = mdata[np.isin(mdata.obs['celltype'], celltypes)].copy()
mdata.obs['celltype'] = mdata.obs['celltype'].cat.remove_unused_categories()

# Downsample
if n_sample > 0:
    n_sample = np.min([n_sample, mdata.obs.shape[0]])
    barcodes = mdata.obs.sample(n=n_sample, random_state=seed, replace=False).index
    mdata = mdata[barcodes, :].copy()

# Extract
rna = mdata.mod['rna']
atac = mdata.mod['atac']

# Make sure enough features
rna = rna[:, np.sum(rna.X.toarray() != 0., axis=0) > 3].copy()
atac = atac[:, np.sum(atac.X.toarray() != 0., axis=0) > 3].copy()

# Normalize
rna.layers['counts'] = rna.X.copy()
atac.layers['counts'] = atac.X.copy()
sc.pp.normalize_total(rna, target_sum=1e4)
sc.pp.log1p(rna)
sc.pp.normalize_total(atac, target_sum=1e4)
sc.pp.log1p(atac)

# HVG
def filter_hvg(adata, n_hvg):
    sc.pp.highly_variable_genes(adata, batch_key='batch')
    hvg = adata.var.sort_values('highly_variable_nbatches', ascending=False).head(n_hvg).index
    del adata.var['highly_variable']
    del adata.var['means']
    del adata.var['dispersions']
    del adata.var['dispersions_norm']
    del adata.var['highly_variable_nbatches']
    del adata.var['highly_variable_intersection']
    del adata.uns
    del adata.obs['batch']
    return hvg.values.astype('U')


rna.obs['batch'] = mdata.obs['batch']
atac.obs['batch'] = mdata.obs['batch']
hvg = filter_hvg(rna, n_hvg)
hvr = filter_hvg(atac, n_hvr)
rna = rna[:, np.isin(rna.var_names.values.astype('U'), hvg)].copy()
atac = atac[:, np.isin(atac.var_names.values.astype('U'), hvr)].copy()

# Filter cells and intersect
rna = rna[(rna.X.A != 0).sum(1) > 3, :].copy()
atac = atac[(atac.X.A != 0).sum(1) > 3, :].copy()
obs_inter = atac.obs_names.intersection(rna.obs_names)
rna = rna[obs_inter].copy()
atac = atac[obs_inter].copy()

# Update mdata
mdata = mdata[obs_inter, :].copy()
mdata.mod['rna'] = rna
mdata.mod['atac'] = atac

# Infer latent space
mdata.obsm['X_spectral'] = snap.tl.multi_spectral([rna, atac], features=None)[1]

# Integrate
n_samples = mdata.obs['batch'].unique().size
if n_samples > 1:
    sce.pp.harmony_integrate(
        mdata,
        key='batch',
        basis='X_spectral',
        adjusted_basis='X_spectral',
        max_iter_harmony=30
    )

# Umap
sc.pp.neighbors(mdata, use_rep="X_spectral")
sc.tl.umap(mdata)

# Clean
del mdata.obsp

# Desparsify
if issparse(rna.X):
    rna.X = rna.X.A
if issparse(atac.X):
    atac.X = atac.X.A

# Update mdata
mdata.mod['rna'] = rna
mdata.mod['atac'] = atac
mdata.update()

if root != 'None':
    sc.pp.neighbors(mdata, use_rep='X_spectral')
    sc.tl.paga(mdata, groups='celltype')
    mdata.uns['iroot'] = np.flatnonzero(mdata.obs['celltype']  == root)[0]
    sc.tl.dpt(mdata)
    mdata.uns = dict()
    del mdata.obsm['X_diffmap']
    del mdata.obsp

# Save
mdata.write(path_output)
