import scanpy as sc
import argparse
import os

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--input', required=True)
parser.add_argument('-o','--output', required=True)
args = vars(parser.parse_args())

inp = args['input']
out = args['output']

# Read object
adata = sc.read_h5ad(inp)

# Split data by modality
gex = adata[:, adata.var['feature_types'] == 'GEX']
atac = adata[:, adata.var['feature_types'] == 'ATAC']

# Write
gex.write(os.path.join(out, 'gex.h5ad'))
gex.write(os.path.join(out, 'atac.h5ad'))
