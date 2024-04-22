# remake loom files for AUC   anakysis of regulons

# load libraries
import os
import numpy as np
import scanpy as sc
import loompy as lp
import muon as mu
import argparse

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--file', type=str)
parser.add_argument('-o', '--output', type=str, default=None)
args = vars(parser.parse_args())

# input h5ad file
input_h5mu = args["file"]
# output loom file
if args["output"] is None:
    output_loom = input_h5mu.replace('.h5mu', 'rna.loom')
else:
    output_loom = args["output"]

# directories and variables
#wdir = '/pasteur/appa/scratch/rtrimbou/3omics/'
#os.chdir(wdir)

### Transformation ####


#read muon file
mudata = mu.read_h5mu(input_h5mu)

# extract data
adata = mudata['rna']
print(adata)

# make loom file
row_attrs = { 
    "Gene": np.array(adata.var.index),
}
col_attrs = { 
    "CellID":  np.array(adata.obs.index),
    "nGene": np.array(np.sum(adata.X.transpose() > 0, axis=0)).flatten(),
    "nUMI": np.array(np.sum(adata.X.transpose(), axis=0)).flatten(),
}

lp.create(output_loom,
          adata.X.transpose(),
          row_attrs,
          col_attrs )
