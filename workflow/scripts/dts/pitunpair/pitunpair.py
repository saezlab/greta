import os
import scanpy as sc
from pathlib import Path
import pandas as pd
import numpy as np
import anndata as ad
import mudata as md
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-c','--path_geneids', required=True)
parser.add_argument('-e','--path_peaks', required=True)
parser.add_argument('-f','--path_output', required=True)
parser.add_argument('-g','--path_expr', required=True)
parser.add_argument('-i', '--path_barmap', required=True)
args = vars(parser.parse_args())

path_barmap = args['path_barmap']
path_geneids = args['path_geneids']
path_peaks = args['path_peaks']
path_output = args['path_output']
path_expr = args['path_expr']

# Read gene ids
geneids = pd.read_csv(path_geneids).set_index('symbol')['id'].to_dict()

# Read barmap
barmap = pd.read_csv(path_barmap, index_col=0)
barmap.index.name = None

# Read data
rna = sc.read_10x_h5(path_expr, genome="GRCh38")
del rna.obs
rna.var.index.name = None

# Filter RNA data based on barmap
rna = rna[barmap['RNA'].values, :]
print(barmap)
rna.obs_names = barmap.index

# Filter faulty gene symbols
ensmbls = np.array([geneids[g] if g in geneids else '' for g in rna.var_names])
msk = ensmbls != ''
rna = rna[:, msk].copy()

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
del rna.var

# Read atac data
atac = ad.read_h5ad(path_peaks)

# Filter ATAC data based on barmap and RNA
atac = atac[rna.obs_names, :]

# Create mdata
mdata = md.MuData(
    {'rna': rna, 'atac': atac,},
    obs=barmap
)

# Write
mdata.write(path_output)
