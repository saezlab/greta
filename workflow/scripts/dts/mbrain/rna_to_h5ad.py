#!/usr/bin/env python3
import sys
from pathlib import Path
import pandas as pd
import scipy.sparse as sp
import anndata as ad

# args
rna_path       = Path(sys.argv[1])
annot_path     = Path(sys.argv[2])
gene_map_path  = Path(sys.argv[3])
celltypes_path = Path(sys.argv[4])
out_path       = Path(sys.argv[5])

# annot (index = ATAC barcodes with sample prefix)
annot = pd.read_csv(annot_path, index_col=0)

# RNA counts
df = pd.read_csv(
    rna_path, sep="\t", index_col=0,
    compression=("gzip" if str(rna_path).endswith(".gz") else None)
)
df.columns = [c.replace(",", ".") for c in df.columns]

# map RNA → ATAC using celltype.txt
ct = pd.read_csv(celltypes_path, sep="\t")  # atac.bc, rna.bc, celltype
ct["rna.bc"]  = ct["rna.bc"].astype(str).str.strip().str.replace(",", ".", regex=False)
ct["atac.bc"] = ct["atac.bc"].astype(str).str.strip().str.replace(",", ".", regex=False)
rna2atac = dict(zip(ct["rna.bc"], ct["atac.bc"]))
df.rename(columns=rna2atac, inplace=True)

# add sample prefix to ATAC barcodes to match annot index
sample_prefix = annot.index[0].split("_")[0] if len(annot) and "_" in annot.index[0] else None
if sample_prefix:
    df.columns = [f"{sample_prefix}_{c}" for c in df.columns]
    print(f"[INFO] Columns renamed to ATAC with prefix '{sample_prefix}'")

# gene map
gene_map = pd.read_csv(gene_map_path, sep="\t")

# build lookups
gene_lookup = {}
gene_lookup_lower = {}
for row in gene_map.itertuples(index=False):
    symbol = str(row.symbol) if pd.notna(row.symbol) else ""
    if not symbol:
        continue
    for ident in [row.ensembl_gene_id, row.symbol]:
        if pd.notna(ident) and str(ident):
            s = str(ident)
            gene_lookup.setdefault(s, symbol)
            gene_lookup_lower.setdefault(s.lower(), symbol)
    if pd.notna(row.synonyms) and row.synonyms:
        for syn in str(row.synonyms).split("|"):
            syn = syn.strip()
            if syn:
                gene_lookup.setdefault(syn, symbol)
                gene_lookup_lower.setdefault(syn.lower(), symbol)

# map gene rows
original_index = df.index.astype(str)
mapped_genes, unmapped = [], []
replaced = 0
for g in original_index:
    if g in gene_lookup:
        new = gene_lookup[g]
        if new != g: replaced += 1
        mapped_genes.append(new)
    elif g.lower() in gene_lookup_lower:
        new = gene_lookup_lower[g.lower()]
        replaced += 1
        mapped_genes.append(new)
    else:
        mapped_genes.append(g)
        unmapped.append(g)
df.index = pd.Index(mapped_genes)

# report unmapped
if unmapped:
    uf = out_path.parent / "unmapped_genes.txt"
    with open(uf, "w") as f:
        f.write(f"# {len(set(unmapped))} genes could not be mapped\n")
        for g in sorted(set(unmapped)):
            f.write(f"{g}\n")
    print(f"WARNING: {len(set(unmapped))} unique genes could not be mapped. See: {uf}")

# collapse duplicates
if df.index.duplicated().any():
    df = df.groupby(level=0, sort=False).sum()

# align to annot (ATAC+prefix)
df = df.loc[:, df.columns.isin(annot.index)]
df = df.loc[:, annot.index]

# AnnData
X = sp.csr_matrix(df.values)
adata = ad.AnnData(X=X.T)
adata.obs = annot.copy()
adata.var_names = df.index.astype(str)
adata.obs_names = df.columns.astype(str)

# write + summary
out_path.parent.mkdir(parents=True, exist_ok=True)
adata.write(out_path)

adata = ad.read_h5ad(out_path)
print("=== AnnData summary ===")
print(adata)
print("\n-- obs (head) --")
print(adata.obs.head())
print("\n-- var (head) --")
print(adata.var.head())
print("\n-- X --")
print(type(adata.X), "shape:", adata.X.shape)
try:
    nnz = (adata.X != 0).sum()
    print("nonzeros:", nnz)
except Exception as e:
    print("nonzeros: n/a;", e)
print("\nkeys: layers", list(adata.layers.keys()), "obsm", list(adata.obsm.keys()),
      "varm", list(adata.varm.keys()), "uns", list(adata.uns.keys()))
print(f"\nReplaced gene identifiers for {replaced} entries using {gene_map_path.name}.")