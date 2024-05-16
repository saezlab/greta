import mudata as mu
import numpy as np
import pandas as pd
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--inp_path', required=True)
parser.add_argument('-o','--out_path', required=True)
args = vars(parser.parse_args())

inp_path = args['inp_path']
out_path = args['out_path']

# Read
mdata = mu.read(inp_path)
atac = mdata.mod['atac'].copy()
rna = mdata.mod['rna'].copy()

# Shuffle genes
rng = np.random.default_rng(seed=42)
genes = rng.choice(rna.var_names, rna.var_names.size, replace=False)
rna.var_names = genes

# Shuffle cres
cres = rng.choice(atac.var_names, atac.var_names.size, replace=False)
atac.var_names = cres

# Shuffle obs
obs = mdata.obs.copy()
for col in obs.columns:
    obs[col] = rng.choice(obs[col], obs.shape[0], replace=False)
mdata.obs = obs

# Update
mdata.mod['rna'] = rna
mdata.mod['atac'] = atac

# Write
mdata.write(out_path)
