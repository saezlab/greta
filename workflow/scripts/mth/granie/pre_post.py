import decoupler as dc
import pandas as pd
import numpy as np
import mudata as mu
import scipy.sparse as ss
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--path_input', required=True)
parser.add_argument('-o','--path_out', required=True)
args = vars(parser.parse_args())

path_input = args['path_input']
path_out = args['path_out']

# Read data
print(path_input)
mdata = mu.read(path_input)

# Remove all equal features
msk = np.any(np.diff(mdata.mod['rna'].X, axis=0), axis=0)
rna = mdata.mod['rna'][:, msk].copy()

msk = np.any(np.diff(mdata.mod['atac'].X, axis=0), axis=0)
atac = mdata.mod['atac'][:, msk].copy()

# Save
obs=mdata.obs.copy()
mdata = mu.MuData({
    'rna': rna,
    'atac': atac,
})
mdata.obs = obs

# Write
mdata.write(path_out)
