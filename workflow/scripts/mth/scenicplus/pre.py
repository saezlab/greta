import mudata as mu
import scipy.sparse as scs
import scanpy as sc
import sys


path_ann = sys.argv[1]
path_scn = sys.argv[2]
path_out = sys.argv[3]

# Read
ann = mu.read(path_ann)
scn = mu.read(path_scn)
scn.var.index = scn.var_names.str.replace(':', '-')

# Match
inter_var = ann.var_names.intersection(scn.var_names)
inter_obs = ann.obs_names.intersection(scn.obs_names)
ann = ann[inter_obs, inter_var].copy()
scn = scn[inter_obs, inter_var]

# Update atac counts with topic ones
ann.mod['atac'].layers['counts'] = scs.csr_matrix(scn.mod['scATAC'].X)
ann.mod['atac'].X = scn.mod['scATAC'].X

# Write
ann.write(path_out)
