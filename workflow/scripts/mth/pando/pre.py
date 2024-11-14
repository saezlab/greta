import pandas as pd
import numpy as np
import mudata as mu
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--path_input', required=True)
parser.add_argument('-p','--path_peaks', required=True)
parser.add_argument('-m','--path_matches', required=True)
parser.add_argument('-o','--path_out', required=True)
args = vars(parser.parse_args())

path_input = args['path_input']
path_peaks = args['path_peaks']
path_matches = args['path_matches']
path_out = args['path_out']

# Read
mdata = mu.read(path_input)
df = pd.read_csv(path_peaks)
matches = pd.read_csv(path_matches).iloc[:, 0].values - 1

# Format peaks
new_peaks = (df['seqnames'].astype(str) + '-' + df['start'].astype(str) + '-' + df['end'].astype(str)).values.astype('U')

# Filter
atac = mdata.mod['atac'][:, matches].copy()
atac.var_names = new_peaks
msk = np.sum(atac.X, axis=1) != 0
atac = atac[msk, :].copy()

# Remove missmatched obs
rna = mdata.mod['rna'].copy()
inter = np.intersect1d(rna.obs_names, atac.obs_names)
x_spectral = mdata[inter, :].obsm['X_spectral'].copy()
x_umap = mdata[inter, :].obsm['X_umap'].copy()
obs = mdata.obs.copy()
obs = obs.loc[inter, :]
mdata = mu.MuData(
    {
    'rna': rna[inter, :].copy(),
    'atac': atac[inter, :].copy(),
    }
)
mdata.obsm['X_spectral'] = x_spectral
mdata.obsm['X_umap'] = x_umap
mdata.obs = obs

# Write
mdata.write(path_out)
