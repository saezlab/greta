import scanpy as sc
import argparse
import os
import numpy as np

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--input', required=True)
parser.add_argument('-p','--plot', required=True)
parser.add_argument('-g','--gex', required=True)
parser.add_argument('-a','--atac', required=True)
args = vars(parser.parse_args())

inp = args['input']
plot = args['plot']
p_gex = args['gex']
p_atac = args['atac']

# Read object
adata = sc.read_h5ad(inp)

# Remove cells that are not in the trajectory
adata = adata[~np.isnan(adata.obs['GEX_pseudotime_order']).values]

# Split data by modality
gex = adata[:, adata.var['feature_types'] == 'GEX']
atac = adata[:, adata.var['feature_types'] == 'ATAC']

# Get log-normalized counts
gex.X = gex.layers['counts']
sc.pp.normalize_total(gex, target_sum=1e4)
sc.pp.log1p(gex)

# Get UMAP
sc.pp.highly_variable_genes(gex)
gex.raw = gex
gex = gex[:, gex.var.highly_variable]
sc.pp.scale(gex)
sc.tl.pca(gex)
sc.pp.neighbors(gex)
sc.tl.umap(gex)
sc.tl.leiden(gex)

# Remove noisy cells
gex = gex[~np.isin(gex.obs['leiden'], ['11', '12'])]
atac = atac[~np.isin(gex.obs['leiden'], ['11', '12'])]

# Plot results
fig, ax = plt.subplots(2, 3, figsize=(12,3), facecolor='white', tight_layout=True)
ax = ax.flatten()
sc.pl.pca(gex, color=['leiden'], ax=ax[0], show=False)
sc.pl.pca(gex, color=['cell_type'], ax=ax[1], show=False)
sc.pl.pca(gex, color=['GEX_pseudotime_order'], ax=ax[2], show=False)
sc.pl.umap(gex, color=['leiden'], ax=ax[3], show=False)
sc.pl.umap(gex, color=['cell_type'], ax=ax[4], show=False)
sc.pl.umap(gex, color=['GEX_pseudotime_order'], ax=ax[5], show=False)
fig.savefig(plot)

# Write
gex.write(p_gex)
atac.write(p_atac)
