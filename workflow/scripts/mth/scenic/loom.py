import os
import numpy as np
import pandas as pd
import mudata as mu
import loompy as lp
import argparse


# Init args
parser = argparse.ArgumentParser()

parser.add_argument('-i','--data', required=True)
parser.add_argument('-o','--path_out', required=True)
args = vars(parser.parse_args())

path_input = args['data']
path_out = args['path_out']

# Extract raw counts data and assign labels
mdata = mu.read(path_input)
adata = mdata.mod['rna'].copy()
adata.layers['lognorm'] = adata.X.copy()
adata.X = adata.layers['counts'].copy()
adata.obs['celltype'] = mdata.obs['celltype']
adata.obsm['X_pca'] = mdata.obsm['X_spectral']

# create basic row and column attributes for the loom file:
row_attrs = {
    "Gene": np.array(adata.var_names) ,
}
col_attrs = {
    "CellID": np.array(adata.obs_names) ,
    "nGene": np.array(np.sum(adata.X.transpose() > 0, axis=0)).flatten() ,
    "nUMI": np.array(np.sum(adata.X.transpose(), axis=0)).flatten() ,
}
lp.create(path_out, adata.X.transpose(), row_attrs, col_attrs)
