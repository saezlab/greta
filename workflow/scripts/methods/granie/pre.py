import decoupler as dc
import pandas as pd
import numpy as np
import mudata as mu
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--path_input', required=True)
parser.add_argument('-g','--path_geneids', required=True)
parser.add_argument('-o','--path_out', required=True)
args = vars(parser.parse_args())

path_input = args['path_input']
path_geneids = args['path_geneids']
path_out = args['path_out']

# Read rna adata
mdata = mu.read(path_input)
rna = mdata.mod['rna'].copy()

# Remove genes that have no ENSEMBL id
geneids = pd.read_csv(path_geneids).set_index('symbol')['id'].to_dict()
ensmbls = np.array([geneids[g] if g in geneids else '' for g in rna.var_names])
msk = ensmbls != ''
rna = rna[:, msk].copy()

# Psbulk rna
rna.obs['celltype'] = mdata.obs['celltype']
rna = dc.get_pseudobulk(
    adata=rna,
    sample_col='celltype',
    groups_col=None,
    layer='counts',
    mode='sum',
    min_cells=10,
    min_counts=100,
)
for col in rna.obs.columns:
    del rna.obs[col]
del rna.layers['psbulk_props']

# Psbulk atac
atac = mdata.mod['atac'].copy()
atac.obs['celltype'] = mdata.obs['celltype']
atac = dc.get_pseudobulk(
    adata=atac,
    sample_col='celltype',
    groups_col=None,
    layer='counts',
    mode='sum',
    min_cells=10,
    min_counts=100,
)
for col in atac.obs.columns:
    del atac.obs[col]
del atac.layers['psbulk_props']

# Intersect and generate new object
inter = np.intersect1d(rna.obs_names, atac.obs_names)
mdata = mu.MuData({
    'rna': rna[inter, :].copy(),
    'atac': atac[inter, :].copy(),
})

# Write
mdata.write(path_out)
