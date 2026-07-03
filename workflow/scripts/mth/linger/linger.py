import argparse
import mudata as mu
import pandas as pd
import time
import os

from LingerGRN.pseudo_bulk import *

def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-g','--linger_GRN', required=True)
parser.add_argument('-d','--out_dir', required=True)
parser.add_argument('-m','--path_mdata', required=True)
parser.add_argument('-v','--version', required=True)
parser.add_argument('-gen','--genome', required=True)
parser.add_argument('-md','--mode', required=True)
args = vars(parser.parse_args())

GRNdir = args['linger_GRN'] + "/"
out_dir = args['out_dir'] + "/"
path_mdata = args['path_mdata']
version = args['version']
genome = args['genome']
mode = args['mode']

if mode == 'parallel':
    from preprocess_fast import *
    import LINGER_tr_fast as LINGER_tr
    import LL_net_fast as LL_net
else:
    from LingerGRN.preprocess import *  
    import LingerGRN.LINGER_tr as LINGER_tr
    import LingerGRN.LL_net as LL_net

log(f"Using {mode} mode")

# Read mdata
mdata = mu.read(path_mdata)
adata_RNA = mdata["rna"].copy()
adata_ATAC = mdata["atac"].copy()

# Restore raw counts (pseudo_bulk expects raw counts matrices)
adata_RNA.X = adata_RNA.layers['counts'].copy()
adata_ATAC.X = adata_ATAC.layers['counts'].copy()
log(f"adata_RNA loaded:\n{adata_RNA}\n")
log(f"adata_ATAC loaded:\n{adata_ATAC}\n")

# Add barcode, gene_ids, label columns
adata_RNA.obs['barcode'] = adata_RNA.obs_names
adata_ATAC.obs['barcode'] = adata_ATAC.obs_names
adata_RNA.var['gene_ids'] = adata_RNA.var_names
adata_ATAC.var['gene_ids'] = adata_ATAC.var_names.str.replace('-', ':', n=1)
adata_RNA.obs['label'] = mdata.obs.loc[adata_RNA.obs_names, 'celltype'].values
adata_ATAC.obs['label'] = mdata.obs.loc[adata_ATAC.obs_names, 'celltype'].values

# Add sample column
bc = adata_RNA.obs['barcode'].values
if len(bc[0].split("-")) == 2:
    adata_RNA.obs['sample'] = [int(s.split("-")[1]) for s in bc]
    adata_ATAC.obs['sample'] = [int(s.split("-")[1]) for s in bc]
else:
    adata_RNA.obs['sample'] = 1
    adata_ATAC.obs['sample'] = 1

adata_RNA.var_names_make_unique()
adata_RNA.var['gene_ids'] = adata_RNA.var_names.values  

# Intersect barcodes
shared = adata_RNA.obs_names.intersection(adata_ATAC.obs_names)
adata_RNA = adata_RNA[shared].copy()
adata_ATAC = adata_ATAC[shared].copy()
log(f"adata_RNA for Linger:\n{adata_RNA}\n")
log(f"adata_ATAC for Linger:\n{adata_ATAC}\n")


# Generate pseudo-bulk
log("Generating pseudo-bulk / metacells...")
samplelist = list(set(adata_ATAC.obs['sample'].values))
TG_pseudobulk = pd.DataFrame([])        # TG x metacells
RE_pseudobulk = pd.DataFrame([])        # RE x metacells
singlepseudobulk = adata_RNA.obs['sample'].nunique() > 10

for tempsample in samplelist:
    adata_RNAtemp = adata_RNA[adata_RNA.obs['sample'] == tempsample]
    adata_ATACtemp = adata_ATAC[adata_ATAC.obs['sample'] == tempsample]
    TG_temp, RE_temp = pseudo_bulk(adata_RNAtemp, adata_ATACtemp, singlepseudobulk)
    TG_pseudobulk = pd.concat([TG_pseudobulk, TG_temp], axis=1)
    RE_pseudobulk = pd.concat([RE_pseudobulk, RE_temp], axis=1)
    RE_pseudobulk[RE_pseudobulk > 100] = 100


# Save pseudobulk and adata matrices
log("Saving pseudobulk data...")
os.makedirs(out_dir + 'data/', exist_ok=True)
adata_ATAC.write(out_dir + 'data/adata_ATAC.h5ad')
adata_RNA.write(out_dir + 'data/adata_RNA.h5ad')
TG_pseudobulk.fillna(0).to_csv(out_dir + 'data/TG_pseudobulk.tsv')
RE_pseudobulk.fillna(0).to_csv(out_dir + 'data/RE_pseudobulk.tsv')
pd.DataFrame(adata_ATAC.var['gene_ids']).to_csv(out_dir + 'data/Peaks.txt', header=None, index=None)


# Train the model
log("Preprocessing and training LINGER model...")


""""
# LINGER assumes following dir strucutre: 
.
├── LINGER_data/
│   └── data_bulk/
│
├── LINGER_output/
│   ├── cell_population_TF_RE_binding.txt
│   ├── cell_population_cis_regulatory.txt
│   ├── cell_population_trans_regulatory.txt
│   └── ...
│
└── data/
    ├── Peaks.txt
    ├── TG_pseudobulk.tsv
    ├── RE_pseudobulk.tsv
    ├── adata_RNA.h5ad
    └── adata_ATAC.h5ad
"""

# LINGER_data/   is dbs/lingerGRN/data_bulk  (GRNdir)
# LINGER_output/ is dts/{org}/{dat}/cases/{case}/runs/linger/  (out_dir)
#   ├── data/    is out_dir + data/ (to keep files together)

# since data/ is not at top level anymore, we need to chdir to out_dir
GRNdir = os.path.abspath(GRNdir) + "/"          
out_dir = os.path.abspath(out_dir) + "/"        
os.chdir(out_dir)

print(GRNdir)
print(out_dir)

preprocess(TG_pseudobulk, RE_pseudobulk, GRNdir, genome, version, out_dir)

activef='ReLU'
LINGER_tr.training(GRNdir, version, out_dir, activef, species='Human')

# Generate regulatory networks
log("Generating cell population GRNs...")
LL_net.TF_RE_binding(GRNdir, adata_RNA, adata_ATAC, genome, version, out_dir)
LL_net.cis_reg(GRNdir, adata_RNA, adata_ATAC, genome, version, out_dir)
LL_net.trans_reg(GRNdir, version, out_dir, genome)
log("GRNs generation done")
