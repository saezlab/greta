"""
Process ATAC peaks from Matrix Market format to h5ad
"""
import sys
import pandas as pd
import numpy as np
import anndata as ad
import scipy.io
import scipy.sparse

peaks_mtx = sys.argv[1]
peak_names = sys.argv[2]
peak_barcodes = sys.argv[3]
annot_csv = sys.argv[4]
output_h5ad = sys.argv[5]

print("Loading peaks counts...")
# Read matrix market file
X = scipy.io.mmread(peaks_mtx).T.tocsr()  # Transpose to cells x peaks

# Read peak names and barcodes
peaks = pd.read_csv(peak_names, header=None, names=['peak'])
barcodes = pd.read_csv(peak_barcodes, header=None, names=['barcode'])

print(f"Matrix shape: {X.shape}")
print(f"Peaks: {len(peaks)}")
print(f"Barcodes: {len(barcodes)}")

# Create AnnData object
adata = ad.AnnData(X=X)
adata.obs_names = barcodes['barcode'].values
adata.var_names = peaks['peak'].values

# Read annotation to filter cells
annot = pd.read_csv(annot_csv, index_col=0)
print(f"Annotation cells: {len(annot)}")

# Filter to cells in annotation
adata = adata[adata.obs_names.isin(annot.index), :].copy()
print(f"After filtering: {adata.shape}")

# Format peak names to match expected format: chr-start-end
# Input format might be chr:start-end or chr-start-end already
new_var_names = []
for p in adata.var_names:
    # Handle both : and - separators
    p = p.replace(':', '-')
    parts = p.split('-')
    if len(parts) == 3:
        seq, start, end = parts
        # Make sure end is properly formatted (0-based)
        end = int(end)
        p = f'{seq}-{start}-{end}'
    new_var_names.append(p)

adata.var_names = new_var_names
print(f"First few peak names: {adata.var_names[:5].tolist()}")

# Make sure we have the right cells in the right order (matching annotation)
adata = adata[annot.index, :].copy()

print(f"Final ATAC data: {adata.shape}")
print(f"Writing to {output_h5ad}")
adata.write(output_h5ad)
print("Done!")
