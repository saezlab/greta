import pandas as pd
import numpy as np
import mudata as mu
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--path_input', required=True)
parser.add_argument('-o','--path_out', required=True)
args = vars(parser.parse_args())

path_input = args['path_input']
path_out = args['path_out']

# Read
mdata = mu.read(path_input)

# Filter
atac = mdata.mod['atac'].copy()

# Remove missmatched obs
rna = mdata.mod['rna'].copy()
inter = np.intersect1d(rna.obs_names, atac.obs_names)
mdata = mu.MuData({
    'rna': rna[inter, :].copy(),
    'atac': atac[inter, :].copy(),
})

# Write
mdata.write(path_out)