#!/usr/bin/env python3
import sys
from pathlib import Path
import pandas as pd

barcodes_path = Path(sys.argv[1])
celltypes_path = Path(sys.argv[2])
batch = sys.argv[3]

# Load input files
barcodes = pd.read_csv(barcodes_path, header=None, names=["rna.bc"], sep="\t")
celltypes = pd.read_csv(celltypes_path, sep="\t")  # columns: atac.bc, rna.bc, celltype

barcodes["rna.bc"] = barcodes["rna.bc"].astype(str).str.strip()
celltypes["rna.bc"] = celltypes["rna.bc"].astype(str).str.strip()
celltypes["atac.bc"] = celltypes["atac.bc"].astype(str).str.strip()

# Merge: match RNA barcodes from barcodes.txt to celltypes.txt
# but keep ATAC barcodes in the output
df = barcodes.merge(
    celltypes[["atac.bc", "rna.bc", "celltype"]],
    on="rna.bc",
    how="left"
)

# Use ATAC barcode as the primary ID (first column)
df["barcode"] = batch + "_" + df["atac.bc"]

# Add batch column
df["batch"] = batch

# Keep required columns only
df = df[["barcode", "batch", "celltype"]]
df.columns = ["", "batch", "celltype"]

# Write output
out_path = barcodes_path.parent / "annot.csv"
df.to_csv(out_path, index=False)

print(f"[DONE] Wrote annotation file: {out_path}")
print(f"Total barcodes: {len(df)} | Matched celltypes: {df['celltype'].notna().sum()}")