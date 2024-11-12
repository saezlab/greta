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

# Read rna adata
mdata = mu.read(path_input)
rna = mdata.mod['rna'].copy()

# Psbulk rna
rna.obs['batch'] = mdata.obs['batch']
rna.obs['celltype'] = mdata.obs['celltype']
rna_b_per_c = (
    rna.obs.reset_index()
    .groupby('celltype', as_index=False)['index']
    .agg(list).set_index('celltype')['index']
    .to_dict()
)
rna = dc.get_pseudobulk(
    adata=rna,
    sample_col='batch',
    groups_col='celltype',
    layer='counts',
    mode='sum',
    min_cells=10,
    min_counts=1000,
)
del rna.obs['psbulk_n_cells']
del rna.obs['psbulk_counts']
del rna.layers['psbulk_props']
rna.layers['counts'] = ss.csr_matrix(rna.X.copy())
rna.uns['rna_b_per_c'] = rna_b_per_c

# Psbulk atac
atac = mdata.mod['atac'].copy()
atac.obs['batch'] = mdata.obs['batch']
atac.obs['celltype'] = mdata.obs['celltype']
atac_b_per_c = (
    atac.obs.reset_index()
    .groupby('celltype', as_index=False)['index']
    .agg(list).set_index('celltype')['index']
    .to_dict()
)
atac = dc.get_pseudobulk(
    adata=atac,
    sample_col='batch',
    groups_col='celltype',
    layer='counts',
    mode='sum',
    min_cells=10,
    min_counts=1000,
)
del atac.obs['psbulk_n_cells']
del atac.obs['psbulk_counts']
del atac.layers['psbulk_props']
atac.layers['counts'] = ss.csr_matrix(atac.X.copy())
atac.uns['atac_b_per_c'] = atac_b_per_c

# Intersect and generate new object
inter = np.intersect1d(rna.obs_names, atac.obs_names)
mdata = mu.MuData({
    'rna': rna[inter, :].copy(),
    'atac': atac[inter, :].copy(),
})
mdata.obs = mdata.mod['rna'].obs.copy()
del mdata.mod['rna'].obs
del mdata.mod['atac'].obs

# Write
mdata.write(path_out)
