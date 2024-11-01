import argparse, os, sys
import pandas as pd
import numpy as np
import mudata as md
import dictys


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-m','--mudata_path', required=True)
parser.add_argument('-t','--tmp_path', required=True)
parser.add_argument('-o','--out_path', required=True)
args = vars(parser.parse_args())

mudata_path = args['mudata_path']
tmp_path = args['tmp_path']
out_path = args['out_path']

# Read
mdata = md.read(mudata_path)

# Process rna
pd.DataFrame(
    np.array(mdata.mod['rna'].layers['counts'].todense()).T,
    columns=mdata.mod['rna'].obs.index,
    index=mdata.mod['rna'].var.index
).to_csv(tmp_path, sep="\t", compression="gzip")

dictys.preproc.qc_reads(tmp_path, tmp_path, 50, 10, 0, 200, 100, 0)
rna_df = pd.read_csv(tmp_path, sep='\t', compression="gzip", index_col=0)
genes, barcodes = rna_df.index.values.astype('U'), rna_df.columns.values.astype('U')
rna = mdata.mod['rna']
rna = rna[barcodes, :][:, genes].copy()
rna.X = rna.layers['counts'].todense().A.copy()

# Process atac
atac = mdata.mod['atac']
atac.X = atac.layers['counts'].todense().A.copy()

# Update
mdata.mod['rna'] = rna
mdata.mod['atac'] = atac
mdata.update()

# Write
mdata.write(out_path)
