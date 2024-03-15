import os
import argparse
import scanpy as sc
import muon as mu
from mudata import MuData
import numpy as np


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--input', required=True)
parser.add_argument('-o','--out', required=True)
args = vars(parser.parse_args())

inp = args['input']
out = args['out']

# Read object
atac = sc.read_h5ad(inp)

# Clean object and restructure
del atac.uns
del atac.obsm
cols = [
    'cell_type', 'batch', 'ATAC_pseudotime_order', 'GEX_pseudotime_order', 'Samplename', 'Site', 'DonorNumber',
    'Modality', 'VendorLot', 'DonorID', 'DonorAge', 'DonorBMI', 'DonorBloodType', 'DonorRace', 'Ethnicity',
    'DonorGender', 'QCMeds', 'DonorSmoker'
]
obs = atac.obs[cols].rename(columns={'cell_type': 'celltype'}).copy()
var = atac.var.copy()
del atac.obs
del atac.var

# Extract counts
atac.X = atac.layers['counts']
del atac.layers['counts']

# Split rna and atac
rna = atac[:, var['feature_types'] == 'GEX'].copy()
atac = atac[:, var['feature_types'] == 'ATAC'].copy()

# Format obs and var
rna.obs['pseudotime_order'] = obs['GEX_pseudotime_order']
atac.obs['pseudotime_order'] = obs['ATAC_pseudotime_order']
obs = obs.drop(columns=['GEX_pseudotime_order', 'ATAC_pseudotime_order'])
rna.var['gene_id'] = var['gene_id']
rna.var_names_make_unique()
atac.var_names_make_unique()

# Remove cells witout trajectory
barcodes = rna.obs['pseudotime_order'].dropna().index
rna = rna[barcodes, :].copy()
barcodes = atac.obs['pseudotime_order'].dropna().index
atac = atac[barcodes, :].copy()

# Basic filtering
sc.pp.filter_cells(rna, min_genes=200)
sc.pp.filter_genes(rna, min_cells=3)
sc.pp.filter_cells(atac, min_genes=2000)
sc.pp.filter_genes(atac, min_cells=3)

# Intersect
obs_inter = atac.obs_names.intersection(rna.obs_names)
rna = rna[obs_inter].copy()
atac = atac[obs_inter].copy()

# Remove qc vars
del rna.obs['n_genes']
del rna.var['n_cells']
del atac.obs['n_genes']
del atac.var['n_cells']

# Init mudata
mdata = MuData({'atac': atac, 'rna': rna})
mdata.obs = obs.loc[obs_inter].copy()
mu.pp.intersect_obs(mdata)
mdata.update()

## Write
mdata.write(out)