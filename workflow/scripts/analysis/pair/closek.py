import mudata as mu
import pandas as pd
import numpy as np
import scanpy as sc
from tqdm import tqdm
from sklearn.neighbors import NearestNeighbors
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-a','--mdata_path', required=True)
parser.add_argument('-b','--barmap_path', required=True)
parser.add_argument('-c','--out_path', required=True)
args = parser.parse_args()


# REad
mdata = mu.read(args.mdata_path)
barmap = pd.read_csv(args.barmap_path)

# Format RNA barmap
barmap.loc[:, 'RNA'] = ['smpl_' + b.replace('-1', '') for b in barmap['RNA']]

# Make sure intersection of all
inter = np.intersect1d(np.intersect1d(barmap['ATAC'], barmap['RNA']), mdata.obs_names)
msk = barmap['ATAC'].isin(inter) & barmap['RNA'].isin(inter)
barmap = barmap.loc[msk, :].reset_index(drop=True)
mdata = mdata[inter, :].copy()

def calculate_k(mdata, omic_a, omic_b, barmap):
    # Compute KNN
    omic = mdata.mod[omic_b]
    omic.obsm['X_spectral'] = mdata.obsm['X_spectral'].copy()
    nn = NearestNeighbors(n_neighbors=omic.shape[0], leaf_size=1, n_jobs=1).fit(omic.obsm['X_spectral'])
    indices = nn.kneighbors(omic.obsm['X_spectral'], return_distance=False)
    df = []
    bar_lst = list(omic.obs_names)
    # Find closest k
    for barcode_a, barcode_b in zip(barmap[omic_a.upper()], barmap[omic_b.upper()]):
        i = bar_lst.index(barcode_a)
        knn = indices[i]
        sorted_barcodes = omic.obs_names.values[knn]
        k = list(sorted_barcodes).index(barcode_b) + 1
        df.append([omic_a, barcode_a, k])
    df = pd.DataFrame(df, columns=['anchor', 'barcode', 'k'])
    return df
    
df = []
for pair in [('atac', 'rna'), ('rna', 'atac')]:
    omic_a, omic_b = pair
    df.append(calculate_k(mdata, omic_a, omic_b, barmap))
df = pd.concat(df)

# Write
df.to_csv(args.out_path, index=False)
