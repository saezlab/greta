import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import celloracle as co
import muon as mu
import scipy
import os
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--path_input', required=True)
parser.add_argument('-k','--knn', required=True)
parser.add_argument('-o','--path_out', required=True)
args = vars(parser.parse_args())

path_input = args['path_input']
k = int(args['knn'])
path_out = args['path_out']

# Read rna adata
mdata = mu.read(path_input)

# Extract raw counts data and assign labels
adata = mdata.mod['rna'].copy()
adata.layers['lognorm'] = adata.X.copy()
adata.X = adata.layers['counts'].copy()
adata.obs['celltype'] = mdata.obs['celltype']
adata.obsm['X_pca'] = mdata.obsm['X_spectral']

# Instantiate Oracle object
oracle = co.Oracle()
oracle.import_anndata_as_raw_count(
    adata=adata,
    cluster_column_name="celltype",
    embedding_name="X_pca"
)

# Compute PCA and select top pcs
oracle.perform_PCA()
n_comps = np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_))>0.002))[0][0]
n_comps = min(n_comps, 50)

# Run imputation
oracle.knn_imputation(
    n_pca_dims=n_comps,
    k=k,
    balanced=True,
    b_sight=k*8,
    b_maxl=k*4,
    n_jobs=os.cpu_count(),
)

# Update object with imputet counts
mdata['rna'].X = oracle.adata.layers['imputed_count']

# Write
mdata.write(path_out)
