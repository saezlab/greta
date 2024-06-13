import os
import scanpy as sc
from pathlib import Path
import pandas as pd
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-t','--path_tmp', required=True)
parser.add_argument('-a','--path_annot', required=True)
args = vars(parser.parse_args())

path_tmp = args['path_tmp']
path_annot = args['path_annot']

# Change default cache dir
if not os.path.exists(path_tmp):
    os.mkdir(path_tmp)
_datasets = datasets()
_datasets.path = Path(path_tmp)

# Download
data = sc.read_10x_h5((path_data), genome="GRCh38", gex_only=True)

## Preprocess to get clustering for annotation
# Basic filtering
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Annotate the group of mitochondrial genes as 'mt'
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Filter cells following standard QC criteria.
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]

# Identify the highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

# Regress and scale the data
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)

# Generate PCA features
sc.tl.pca(adata, svd_solver='arpack')

# Restore X to be norm counts
dc.swap_layer(adata, 'log_norm', X_layer_key=None, inplace=True)

# Normalize the data
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers['log_norm'] = adata.X.copy()

# Compute distances in the PCA space, and find cell neighbors
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# Generate UMAP features
sc.tl.umap(adata)

# Run leiden clustering algorithm
sc.tl.leiden(adata, resolution=0.3)


# annotate cell types 
annotation_dict = {
 '0': 'Gonadotropes',
 '1': 'Stem cells',
 '2': 'Pituicytes',
 '3': 'Somatotropes',
 '4': 'Lactotropes',
 '5': 'Thyrotropes',
 '6': 'Corticotropes',
 '7': 'Endothelial cells',
 '8': 'Immune cells',
 '9': 'Pericytes',
 '10': 'Pituicytes',}

adata.obs['celltype'] = [annotation_dict[clust] for clust in adata.obs['leiden']]


# Extract annot
data.obs['batch'] = 'smpl'
data.obs['celltype'] = 'celltype'
annot = data.obs[['batch', 'celltype']]
annot.index.name = None

# Write
annot.to_csv(path_annot)