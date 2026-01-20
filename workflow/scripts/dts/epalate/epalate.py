import mudata as mu
import pandas as pd
import numpy as np
import scanpy as sc
import mudata as mu
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-a','--path_gex', required=True)
parser.add_argument('-b','--path_peaks', required=True)
parser.add_argument('-c','--path_annot', required=True)
parser.add_argument('-d','--path_output', required=True)
args = vars(parser.parse_args())

path_gex = args['path_gex']
path_peaks = args['path_peaks']
path_annot = args['path_annot']
path_output = args['path_output']

# Read
rna = sc.read_h5ad(path_gex)
atac = sc.read_h5ad(path_peaks)
obs = pd.read_csv(path_annot, index_col=0)

# Remove all metadata from RNA (keep only index)
rna.obs = pd.DataFrame(index=rna.obs_names)

# Match
common_cells = rna.obs_names[
    rna.obs_names.isin(atac.obs_names) & rna.obs_names.isin(obs.index)
]
rna = rna[common_cells].copy()
atac = atac[common_cells].copy()
obs = obs.loc[common_cells].copy()
del rna.var

# Create mdata (global obs only)
mdata = mu.MuData(
    {'rna': rna, 'atac': atac},
    obs=obs
)

# Write
mdata.write(path_output)
