#!/usr/bin/env python3
"""
Process downloaded ATAC peaks and counts for lateanagen dataset.
Converts BED peaks and count matrix into h5ad format matching callpeaks.py output.
"""
import sys
from pathlib import Path
import pandas as pd
import scipy.sparse as sp
import anndata as ad

peaks_bed_path = Path(sys.argv[1])  # peaks BED file
counts_path = Path(sys.argv[2])     # counts matrix
annot_path = Path(sys.argv[3])      # annotation CSV
sample_id = sys.argv[4]             # sample name (e.g., 'GSM4156597')
out_path = Path(sys.argv[5])        # output h5ad path

print(f"[INFO] Processing ATAC peaks for {sample_id}...")

# Read annotation
annot = pd.read_csv(annot_path, index_col=0)
print(f"[INFO] Loaded {len(annot)} cells from annotation")

# Read peaks BED file
peaks_df = pd.read_csv(peaks_bed_path, sep="\t", header=None, 
                       names=["chr", "start", "end"], compression="gzip")
print(f"[INFO] Loaded {len(peaks_df)} peaks from BED file")

# Create peak names in format: chr-start-end (matching callpeaks.py output)
# Note: callpeaks.py subtracts 1 from end position, so we do the same
peak_names = [f"{row.chr}-{row.start}-{row.end-1}" for row in peaks_df.itertuples(index=False)]

# Read counts matrix (peaks x cells or cells x peaks - need to check)
counts_df = pd.read_csv(counts_path, sep="\t", compression="gzip", index_col=0)
print(f"[INFO] Loaded counts matrix (raw): {counts_df.shape[0]} rows x {counts_df.shape[1]} columns")

# Check if matrix needs to be transposed
# If we have way more rows than columns, and columns == 0, the file might be cells x peaks (transposed)
# Or if the number of rows matches expected peaks, it's peaks x cells
if counts_df.shape[1] == 0 or counts_df.shape[0] > 10 * len(peaks_df):
    print("[WARNING] Counts matrix appears transposed (cells x peaks). Transposing...")
    counts_df = counts_df.T
    print(f"[INFO] After transpose: {counts_df.shape[0]} rows x {counts_df.shape[1]} columns")

# Convert barcode column names from commas to periods (same as RNA processing)
counts_df.columns = [bc.replace(",", ".") for bc in counts_df.columns]
print(f"[INFO] Processed counts matrix: {counts_df.shape[0]} peaks x {counts_df.shape[1]} cells")

# Prepend sample name to barcodes to match annotation format
# The annotation barcodes are formatted as "sampleID_barcode"
counts_df.columns = [f"{sample_id}_{bc}" if not bc.startswith(sample_id) else bc 
                     for bc in counts_df.columns]

# Filter counts to only include annotated cells
common_cells = annot.index.intersection(counts_df.columns)
print(f"[INFO] Found {len(common_cells)} cells in both annotation and counts")

if len(common_cells) == 0:
    print("[ERROR] No common cells found between annotation and counts!")
    print(f"[DEBUG] Sample annotation index (first 5): {list(annot.index[:5])}")
    print(f"[DEBUG] Sample counts columns (first 5): {list(counts_df.columns[:5])}")
    raise ValueError("No overlapping cells between annotation and counts")

counts_df = counts_df.loc[:, common_cells]
annot = annot.loc[common_cells, :]

# Verify we have the right number of peaks
if len(counts_df) != len(peak_names):
    raise ValueError(f"Mismatch: {len(counts_df)} rows in counts but {len(peak_names)} peaks in BED")

# Reorder to match annotation order
counts_df = counts_df.loc[:, annot.index]

# Create AnnData object (cells x peaks, matching callpeaks.py output)
X = sp.csr_matrix(counts_df.values.T)
adata = ad.AnnData(X=X)
adata.obs_names = annot.index
adata.var_names = peak_names

# Clean up obs and var DataFrames to match callpeaks.py output
# callpeaks.py removes obs and var before saving
del adata.obs
del adata.var

out_path.parent.mkdir(parents=True, exist_ok=True)
adata.write(out_path)

print(f"[INFO] Wrote peaks h5ad to {out_path}")
print(f"[INFO] Shape: {adata.shape[0]} cells x {adata.shape[1]} peaks")
print("[DONE]")
