import os
import pandas as pd
import numpy as np
import mudata as mu
import sys
import scanpy as sc


# Read
path_data = sys.argv[1]
path_out = sys.argv[2]

# Read
mdata = mu.read(path_data)
n_hvg = 4096
rna = mdata.mod['rna']
rna.obs = mdata.obs
sc.pp.highly_variable_genes(rna, batch_key='batch')
hvg = rna.var.sort_values('highly_variable_nbatches', ascending=False).head(n_hvg).index
hvg = pd.DataFrame(hvg.values.reshape(-1, 1), columns=['gene'])
hvg.to_csv(path_out, index=False, header=False, compression="gzip")
