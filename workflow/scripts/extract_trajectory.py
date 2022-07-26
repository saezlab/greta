import pandas as pd
import numpy as np
import muon as mu
from muon import atac as ac
import scanpy as sc
import scanpy.external as sce
import os
import matplotlib.pyplot as plt
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--path_input', required=True)
parser.add_argument('-c','--celltypes', required=True)
parser.add_argument('-p','--plot', required=True)
parser.add_argument('-o','--path_output', required=True)
args = vars(parser.parse_args())

path_input = args['path_input']
celltypes = args['celltypes']
plot = args['plot']
path_output = args['path_output']

# Read object
mdata = mu.read(path_input)

# Process celltypes
celltypes = celltypes.split(';')

# Filter
mdata = mdata[np.isin(mdata.obs['celltype'], celltypes)]

# RNA

## Extract
rna = mdata.mod['rna']
del rna.X
rna.X = rna.layers['counts']

## Normalize
sc.pp.normalize_total(rna, target_sum=1e4)
sc.pp.log1p(rna)

## HVG
sc.pp.highly_variable_genes(rna, min_mean=0.02, max_mean=4, min_disp=0.5)

## PCA
rna.raw = rna
sc.pp.scale(rna, max_value=10)
sc.tl.pca(rna, svd_solver='arpack', use_highly_variable=True)

## Integrate
rna.obs['batch'] = mdata.obs['batch']
sce.pp.harmony_integrate(rna, 'batch', adjusted_basis='X_pca', max_iter_harmony=30)

## UMAP
sc.pp.neighbors(rna, n_neighbors=10, n_pcs=20)
sc.tl.leiden(rna, resolution=.5)
sc.tl.umap(rna, spread=1., min_dist=.5, random_state=11)

# ATAC

## Extract
atac = mdata.mod['atac']
del atac.X
atac.X = atac.layers['counts']

## Normalize
ac.pp.tfidf(atac, scale_factor=1e4)
sc.pp.normalize_per_cell(atac, counts_per_cell_after=1e4)
sc.pp.log1p(atac)

## HVP
sc.pp.highly_variable_genes(atac, min_mean=0.05, max_mean=1.5, min_disp=.5)

## LSI
atac.raw = atac
ac.tl.lsi(atac)
atac.obsm['X_lsi'] = atac.obsm['X_lsi'][:,1:]
atac.varm["LSI"] = atac.varm["LSI"][:,1:]
atac.uns["lsi"]["stdev"] = atac.uns["lsi"]["stdev"][1:]

## Integrate
atac.obs['batch'] = mdata.obs['batch']
sce.pp.harmony_integrate(atac, 'batch', basis='X_lsi', adjusted_basis='X_lsi', max_iter_harmony=30)

## UMAP
sc.pp.neighbors(atac, use_rep="X_lsi", n_neighbors=10, n_pcs=30)
sc.tl.leiden(atac, resolution=.5)
sc.tl.umap(atac, spread=1.5, min_dist=.5, random_state=20)

# Integration

## Intersect
mu.pp.intersect_obs(mdata)
mdata.update()

## Figure
rna.obs['celltype'] = mdata.obs['celltype']
atac.obs['celltype'] = mdata.obs['celltype']
fig, axes = plt.subplots(2,3, figsize=(12,8), dpi=100, facecolor='white', tight_layout=True)
axes = axes.flatten()
sc.pl.umap(rna, color='celltype', ax=axes[0], legend_loc='on data', return_fig=False, show=False)
sc.pl.umap(rna, color='leiden', ax=axes[3], legend_loc='on data', return_fig=False, show=False)
sc.pl.umap(atac, color='celltype', ax=axes[1], legend_loc='on data', return_fig=False, show=False)
sc.pl.umap(atac, color='leiden', ax=axes[4], legend_loc='on data', return_fig=False, show=False)
#sc.pl.umap(mdata, color='celltype', ax=axes[2], legend_loc='on data', return_fig=False, show=False)
#sc.pl.umap(mdata, color='leiden_joint', ax=axes[5], legend_loc='on data', return_fig=False, show=False)
fig.savefig(plot)

# Save
mdata.write(path_output)

