import h5py
import numpy as np
from scipy.sparse import csr_matrix, vstack
from tqdm import tqdm
import mudata as mu
import pandas as pd
import pyranges as pr
import anndata as ad
import sys


def get_pr(index):
    df = index.str.replace(':', '-') 
    df = df.str.split('-').tolist()
    df = pd.DataFrame(df, columns=['Chromosome', 'Start', 'End'])
    return pr.PyRanges(df)


def get_vars(df):
    chrm = df.df['Chromosome'].astype(str)
    strt = df.df['Start'].astype(str)
    end = df.df['End'].astype(str)
    return pd.Index(chrm + ':' + strt + '-' + end)


path_ann = sys.argv[1]
path_scn = sys.argv[2]
path_out = sys.argv[3]

# Read datasets
ann = mu.read(path_ann)
pr_ann = get_pr(ann.mod['atac'].var_names)
with h5py.File(path_scn, 'r') as f:
    # Read names
    scn_obs_names = f['mod']['scATAC']['obs']['_index'][:].astype('U')
    scn_rna_var_names = f['mod']['scRNA']['var']['_index'][:].astype('U')
    scn_atac_var_names = f['mod']['scATAC']['var']['_index'][:].astype('U')

    # Find common obs and vars
    obs_names = ann.obs_names.intersection(scn_obs_names)
    rna_var_names = ann.mod['rna'].var_names.intersection(scn_rna_var_names)
    pr_scn = get_pr(pd.Index(scn_atac_var_names))
    atac_var_names = get_vars(pr_scn.overlap(pr_ann))

    # Map common obs and vars
    obs_map_dict = {v:i for i, v in enumerate(scn_obs_names)}
    obs_msk = np.array([obs_map_dict[v] for v in obs_names])
    var_map_dict = {v:i for i, v in enumerate(scn_rna_var_names)}
    rna_var_msk = np.array([var_map_dict[v] for v in rna_var_names])
    var_map_dict = {v:i for i, v in enumerate(scn_atac_var_names)}
    atac_var_msk = np.array([var_map_dict[v] for v in atac_var_names])
    
    # Read counts
    scn_atac_X = csr_matrix(f['mod']['scATAC']['X'][:, atac_var_msk][obs_msk, :])
    scn_rna_X = f['mod']['scRNA']['X'][:, rna_var_msk][obs_msk, :]


# Build mdata
rna = ad.AnnData(
    X=scn_rna_X,
    obs=pd.DataFrame(index=scn_obs_names[obs_msk]),
    var=pd.DataFrame(index=scn_rna_var_names[rna_var_msk]),
    layers={'counts': ann.mod['rna'][obs_names, scn_rna_var_names].layers['counts']}
)
atac = ad.AnnData(
    X=scn_atac_X,
    obs=pd.DataFrame(index=scn_obs_names[obs_msk]),
    var=pd.DataFrame(index=pd.Index(scn_atac_var_names[atac_var_msk]).str.replace(':', '-')),
    layers={'counts': scn_atac_X}
)
assert (ann.mod['rna'][:, rna_var_names].var_names == rna.var_names).all()
assert (ann.mod['rna'][obs_names[0], :].X == rna[obs_names[0], :].X).all()
mdata = mu.MuData(
    {'rna': rna, 'atac': atac}
)
mdata.obs = ann[obs_names, :].obs
mdata.obsm = ann[obs_names, :].obsm

# Write
mdata.write(path_out)
