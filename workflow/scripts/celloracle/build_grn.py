import pandas as pd
import numpy as np
import celloracle as co
import mudata
import os
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-m','--path_mdata', required=True)
parser.add_argument('-b','--base_grn', required=True)
parser.add_argument('-l','--path_links', required=True)
args = vars(parser.parse_args())

path_mdata = args['path_mdata']
path_base_grn = args['base_grn']
path_links = args['path_links']

# Read rna adata
adata = mudata.read(os.path.join(path_mdata, 'rna'))

# Load base grn
base_GRN = pd.read(path_base_grn)

# Instantiate Oracle object
oracle = co.Oracle()

# Set to raw counts
adata.X = adata.layers['counts']

# Filter by HVG
adata = adata[:, adata.var.highly_variable]

# Instantiate Oracle object
oracle.import_anndata_as_raw_count(adata=adata,
                                   cluster_column_name="celltype",
                                   embedding_name="X_umap")

# You can load TF info dataframe with the following code
oracle.import_TF_data(TF_info_matrix=base_GRN)

# Perform PCA
oracle.perform_PCA()

# KNN imputation
n_cell = oracle.adata.shape[0]
print("cell number is: {0}".format(n_cell))

k = int(0.025 * n_cell)
print('Auto-selected k is: {1}'.format(k))

oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8,
                      b_maxl=k*4, n_jobs=os.cpu_count())

# Compute GRN
links = oracle.get_links(cluster_name_for_GRN_unit="celltype", alpha=10,
                         verbose_level=10)

# Save Links object
links.to_hdf5(file_path=path_links)

