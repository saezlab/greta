"""
Process RNA counts from Matrix Market format to h5ad
"""
import sys
import pandas as pd
import numpy as np
import anndata as ad
import scipy.io
import scipy.sparse
import scanpy as sc

rna_mtx = sys.argv[1]
rna_genes = sys.argv[2]
rna_barcodes = sys.argv[3]
annot_csv = sys.argv[4]
gene_ids = sys.argv[5]
output_h5ad = sys.argv[6]

print("Loading RNA counts...")
# Read matrix market file
X = scipy.io.mmread(rna_mtx).T.tocsr()  # Transpose to cells x genes

# Read genes and barcodes
genes = pd.read_csv(rna_genes, header=None, names=['gene'])
barcodes = pd.read_csv(rna_barcodes, header=None, names=['barcode'])

print(f"Matrix shape: {X.shape}")
print(f"Genes: {len(genes)}")
print(f"Barcodes: {len(barcodes)}")

# Create AnnData object
adata = ad.AnnData(X=X)
adata.obs_names = barcodes['barcode'].values
adata.var_names = genes['gene'].values

# Read annotation to filter cells
annot = pd.read_csv(annot_csv, index_col=0)
print(f"Annotation cells: {len(annot)}")

# Filter to cells in annotation
adata = adata[adata.obs_names.isin(annot.index), :].copy()
print(f"After filtering: {adata.shape}")

# Load gene ID mapping (columns are lowercase: id, symbol)
geneids = pd.read_csv(gene_ids)
print(f"Gene ID file columns: {geneids.columns.tolist()}")
geneids = geneids.set_index('symbol')['id'].to_dict()

# Map genes to Ensembl IDs
ensmbls = np.array([geneids.get(g, '') for g in adata.var_names])
msk = ensmbls != ''
print(f"Genes with Ensembl IDs: {msk.sum()} / {len(adata.var_names)}")

# Filter to genes with IDs
adata = adata[:, msk].copy()
adata.var['gene_ids'] = ensmbls[msk]

# Remove duplicated genes (keep first occurrence)
duplicated_genes = adata.var_names[adata.var_names.duplicated()].unique()
if len(duplicated_genes) > 0:
    print(f"Removing {len(duplicated_genes)} duplicated genes")
    adata = adata[:, ~adata.var_names.isin(duplicated_genes)].copy()

# Make sure we have the right cells in the right order (matching annotation)
adata = adata[annot.index, :].copy()

print(f"Final RNA data: {adata.shape}")
print(f"Writing to {output_h5ad}")
adata.write(output_h5ad)
print("Done!")
