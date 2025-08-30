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
parser.add_argument('-a','--path_gex', required=True)
parser.add_argument('-b','--path_peaks', required=True)
parser.add_argument('-c','--path_annot', required=True)
parser.add_argument('-f','--path_output', required=True)
args = vars(parser.parse_args())

path_gex = args['path_gex']
path_peaks = args['path_peaks']
path_annot = args['path_annot']
path_output = args['path_output']

# Read annots
rna = ad.read_h5ad(path_gex)
atac = ad.read_h5ad(path_peaks)
obs = pd.read_csv(path_annot, index_col=0)

# Create mdata
mdata = md.MuData(
    {'rna': rna, 'atac': atac,},
    obs=obs
)

# Write
mdata.write(path_output)
