import pandas as pd
import numpy as np
import muon as mu
from muon import atac as ac
import scanpy as sc
import scanpy.external as sce
import os
import matplotlib.pyplot as plt
from scipy.sparse import issparse
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--path_input', required=True)
parser.add_argument('-c','--celltypes', required=True)
parser.add_argument('-g','--n_hvg', required=True)
parser.add_argument('-r','--n_hvr', required=True)
parser.add_argument('-o','--path_output', required=True)
args = vars(parser.parse_args())

path_input = args['path_input']
celltypes = args['celltypes']
n_hvg = int(args['n_hvg'])
n_hvr = int(args['n_hvr'])
path_output = args['path_output']

# Read object
mdata = mu.read(path_input)

# Process celltypes
celltypes = celltypes.split(';')

# Filter
mdata = mdata[np.isin(mdata.obs['celltype'], celltypes)].copy()
mdata.obs['celltype'] = mdata.obs['celltype'].cat.remove_unused_categories()

# Extract
rna = mdata.mod['rna']
atac = mdata.mod['atac']

# Normalize
rna.layers['counts'] = rna.X.copy()
atac.layers['counts'] = atac.X.copy()
sc.pp.normalize_total(rna, target_sum=1e4)
sc.pp.log1p(rna)
ac.pp.tfidf(atac, scale_factor=1e4)

# HVG
rna.obs['batch'] = mdata.obs['batch']
atac.obs['batch'] = mdata.obs['batch']

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
    adata = adata[:, hvg].copy()
    return adata

rna = filter_hvg(rna, n_hvg=n_hvg)
atac = filter_hvg(atac, n_hvg=n_hvr)

# Filter cells and intersect
rna = rna[(rna.X.A != 0).sum(1) > 3, :].copy()
atac = atac[(atac.X.A != 0).sum(1) > 3, :].copy()
obs_inter = atac.obs_names.intersection(rna.obs_names)
rna = rna[obs_inter].copy()
atac = atac[obs_inter].copy()

# TODO: add downsampling
# Subset if needed
#if n_downsample > 0:
#    index = rna.obs.sample(frac=1, random_state=1, replace=False).index
#    barcodes = mdata.obs.groupby(['celltype', 'batch']).head(n_downsample).index
#    rna = rna[barcodes, :].copy()
#    atac = atac[barcodes, :].copy()


## PCA
rna.layers['norm'] = rna.X.copy()
sc.pp.scale(rna, max_value=10)
sc.tl.pca(rna, svd_solver='arpack')
rna.X = rna.layers['norm'].copy()
del rna.layers['norm']

## Integrate
sce.pp.harmony_integrate(rna, 'batch', adjusted_basis='X_pca', max_iter_harmony=30)

## UMAP
sc.pp.neighbors(rna)
sc.tl.umap(rna)

## LSI
atac.layers['norm'] = atac.X.copy()
ac.tl.lsi(atac)
atac.obsm['X_lsi'] = atac.obsm['X_lsi'][:, 1:]
atac.varm["LSI"] = atac.varm["LSI"][:, 1:]
atac.X = atac.layers['norm'].copy()
del atac.layers['norm']
del atac.uns

## Integrate
sce.pp.harmony_integrate(atac, 'batch', basis='X_lsi', adjusted_basis='X_lsi', max_iter_harmony=30)

## UMAP
sc.pp.neighbors(atac, use_rep="X_lsi")
sc.tl.umap(atac)

# Clean
del rna.obs['batch']
del rna.var['mean']
del rna.var['std']
del rna.uns
del rna.varm
del rna.obsp
del atac.obs['batch']
del atac.uns
del atac.varm
del atac.obsp

# Desparsify
if issparse(rna.X):
    rna.X = rna.X.A
if issparse(atac.X):
    atac.X = atac.X.A

# Update mdata
mdata.mod['rna'] = rna
mdata.mod['atac'] = atac
mdata.update()

# Save
mdata.write(path_output)
