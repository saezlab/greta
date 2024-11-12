import os
import scanpy as sc
import snapatac2 as snap
from snapatac2.datasets import _datasets, datasets
from pathlib import Path
import pandas as pd
import numpy as np
import anndata as ad
import mudata as md
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-a','--path_gex', nargs='+', required=True)
parser.add_argument('-b','--path_peaks', required=True)
parser.add_argument('-c','--path_annot', required=True)
parser.add_argument('-d','--path_geneids', required=True)
parser.add_argument('-f','--path_output', required=True)
args = vars(parser.parse_args())

path_gex = args['path_gex']
path_peaks = args['path_peaks']
path_annot = args['path_annot']
path_geneids = args['path_geneids']
path_output = args['path_output']

print(path_gex)

# Read annots
obs = pd.read_csv(path_annot, index_col=0)
obs.index = [sample_id + '_' + barcode for barcode, sample_id in zip(obs.index, obs['batch'])]

# Read individual h5 files into list of adatas
rnas = [sc.read_10x_h5(p) for p in path_gex]
# Reformat obs names of adatas
for r in rnas:
    r.var_names_make_unique()
    r.obs_names_make_unique()
    r.obs_names = [b.split('-1')[0] for b in r.obs_names]
    
sample_ids = [os.path.basename(p).split('_')[0]for p in path_gex]

for i in range(0, len(sample_ids)):
    rnas[i].obs_names = sample_ids[i] + '_' + rnas[i].obs_names

# Concatenate list of adatas into a single object
rna = ad.concat(rnas)
rna.var['gene_ids'] = rna.var_names


# Read gene ids
geneids = pd.read_csv(path_geneids).set_index('symbol')['id'].to_dict()

# Filter faulty gene symbols
ensmbls = np.array([geneids[g] if g in geneids else '' for g in rna.var_names])
msk = ensmbls != ''
rna = rna[:, msk].copy()
rna.var_names = ensmbls[msk]

# Basic QC
sc.pp.filter_cells(rna, min_genes=100)
sc.pp.filter_genes(rna, min_cells=3)
del rna.obs['n_genes']

# Remove duplicated genes based on num of cells
to_remove = []
for dup in rna.var.index[rna.var.index.duplicated()]:
    tmp = rna.var.loc[dup]
    max_idx = tmp.set_index('gene_ids')['n_cells'].idxmax()
    to_remove.extend(tmp['gene_ids'][tmp['gene_ids'] != max_idx].values)
rna = rna[:, ~rna.var['gene_ids'].isin(to_remove)].copy()
del rna.uns
del rna.obsm
del rna.obsp
del rna.obs
del rna.var

# Read atac data
atac = ad.read_h5ad(path_peaks)
rna = rna[atac.obs_names].copy()
atac = atac[rna.obs_names].copy()

# Create mdata
mdata = md.MuData(
    {'rna': rna, 'atac': atac,},
    obs=obs
)

# Write
mdata.write(path_output)