"""
Combine RNA and ATAC h5ad files into MuData for retina dataset
"""
import mudata as mu
import pandas as pd
import scanpy as sc
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-a', '--path_gex', required=True)
parser.add_argument('-b', '--path_peaks', required=True)
parser.add_argument('-c', '--path_annot', required=True)
parser.add_argument('-d', '--path_output', required=True)
args = vars(parser.parse_args())

path_gex = args['path_gex']
path_peaks = args['path_peaks']
path_annot = args['path_annot']
path_output = args['path_output']

# Read
rna = sc.read_h5ad(path_gex)
atac = sc.read_h5ad(path_peaks)
obs = pd.read_csv(path_annot, index_col=0)
# Match
atac = atac[rna.obs_names].copy()
# Convert matrices to match skin dataset format
rna.X = rna.X.tocsr().astype('int32')
atac.X = atac.X.tocsr().astype('uint32')
# Remove gene_ids column from RNA .var
if 'gene_ids' in rna.var.columns:
    rna.var = rna.var.drop(columns=['gene_ids'])
# Remove obs columns from modalities
rna.obs = rna.obs[[]]
atac.obs = atac.obs[[]]
# Create mdata
mdata = mu.MuData(
    {'rna': rna, 'atac': atac,},
    obs=obs
)
# Write
mdata.write(path_output)
