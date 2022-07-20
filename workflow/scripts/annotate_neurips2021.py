import scanpy as sc
import scanpy.external as sce
import muon as mu
from muon import atac as ac
import argparse
import os
from mudata import MuData
import numpy as np
import matplotlib.pyplot as plt


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--input', required=True)
parser.add_argument('-p','--plot', required=True)
parser.add_argument('-g','--use_gpu', required=True)
parser.add_argument('-o','--out', required=True)
args = vars(parser.parse_args())

inp = args['input']
plot = args['plot']
use_gpu = args['use_gpu']
if use_gpu == 'True':
    use_gpu = True
else:
    use_gpu = False
out = args['out']

# Read object
mdata = sc.read_h5ad(inp)

# Clean object and restructure
del mdata.uns
del mdata.obsm
mdata = MuData({'atac': mdata[:, mdata.var['feature_types'] == 'ATAC'],
                'rna': mdata[:, mdata.var['feature_types'] == 'GEX']})
rna = mdata.mod['rna']
atac = mdata.mod['atac']
for col in rna.obs.columns:
    if not col.startswith('GEX_') and not col.startswith('ATAC_'):
        mdata.obs[col] = rna.obs[col]
mdata.obs.rename(columns = {'cell_type': 'celltype'}, inplace = True)

def filter_names(data, prefix):
    for col in data.obs.columns:
        if not col.startswith(prefix):
            del data.obs[col]
        else:
            data.obs.rename(columns = {col: col.split(prefix)[1]}, inplace = True)

filter_names(rna, 'GEX_')
filter_names(atac, 'ATAC_')
del rna.X
rna.X = rna.layers['counts']
del atac.X
atac.X = atac.layers['counts']

# RNA

## QC
rna.var['mt'] = rna.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(rna, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
mu.pp.filter_var(rna, 'n_cells_by_counts', lambda x: x >= 3)
mu.pp.filter_obs(rna, 'n_genes_by_counts', lambda x: (x >= 200) & (x < 5000))
mu.pp.filter_obs(rna, 'total_counts', lambda x: x < 15000)
mu.pp.filter_obs(rna, 'pct_counts_mt', lambda x: x < 15)

## Normalize
sc.pp.normalize_total(rna, target_sum=1e4)
sc.pp.log1p(rna)

## HVG
sc.pp.highly_variable_genes(rna, min_mean=0.02, max_mean=4, min_disp=0.5)

## PCA
rna.raw = rna
rna = rna[:, rna.var['highly_variable']]
sc.pp.scale(rna, max_value=10)
sc.tl.pca(rna, svd_solver='arpack')

## Integrate
rna.obs['batch'] = mdata.obs['batch']
sce.pp.harmony_integrate(rna, 'batch', adjusted_basis='X_pca', max_iter_harmony=30)

## UMAP
sc.pp.neighbors(rna, n_neighbors=10, n_pcs=20)
sc.tl.leiden(rna, resolution=.5)
sc.tl.umap(rna, spread=1., min_dist=.5, random_state=11)

# ATAC

## QC
sc.pp.calculate_qc_metrics(atac, percent_top=None, log1p=False, inplace=True)
mu.pp.filter_var(atac, 'n_cells_by_counts', lambda x: x >= 10)
mu.pp.filter_obs(atac, 'n_genes_by_counts', lambda x: (x >= 2000) & (x <= 15000))
mu.pp.filter_obs(atac, 'total_counts', lambda x: (x >= 4000) & (x <= 40000))

## Normalize
ac.pp.tfidf(atac, scale_factor=1e4)
sc.pp.normalize_per_cell(atac, counts_per_cell_after=1e4)
sc.pp.log1p(atac)

## HVP
sc.pp.highly_variable_genes(atac, min_mean=0.05, max_mean=1.5, min_disp=.5)

## LSI
atac.raw = atac
atac = atac[:, atac.var['highly_variable']]
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

## MOFA
#mu.tl.mofa(mdata, groups_label='batch', verbose=True, use_raw=True)

## UMAP
#sc.pp.neighbors(mdata, use_rep="X_mofa")
#sc.tl.leiden(mdata, key_added='leiden_joint', resolution=0.5)
#sc.tl.umap(mdata, min_dist=.2, spread=1., random_state=10)

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

## Write
mdata.write(out)
