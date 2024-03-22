import pandas as pd
import numpy as np
import mudata as mu
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--path_input', required=True)
parser.add_argument('-p','--path_peaks', required=True)
parser.add_argument('-o','--path_out', required=True)
args = vars(parser.parse_args())

path_input = args['path_input']
path_peaks = args['path_peaks']
path_out = args['path_out']

# Read
mdata = mu.read(path_input)
df = pd.read_csv(path_peaks)

# Format peaks
peaks = df['seqnames'].astype(str) + '-' + df['start'].astype(str) + '-' + df['end'].astype(str)

# Filter
atac = mdata.mod['atac'][:, peaks].copy()
msk = np.sum(atac.X, axis=1) != 0
atac = atac[msk, :].copy()

# Remove missmatched obs
rna = mdata.mod['rna'].copy()
inter = np.intersect1d(rna.obs_names, atac.obs_names)
mdata = mu.MuData({
    'rna': rna[inter, :].copy(),
    'atac': atac[inter, :].copy(),
})

# Write
mdata.write(path_out)
