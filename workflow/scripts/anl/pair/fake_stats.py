import mudata as mu
import pandas as pd
import numpy as np
import scanpy as sc
from tqdm import tqdm
import scipy.stats as st
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import adjusted_rand_score as ari
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-a','--mdata_path', required=True)
parser.add_argument('-b','--barmap_path', required=True)
parser.add_argument('-c','--knn_path', required=True)
parser.add_argument('-d','--cor_path', required=True)
parser.add_argument('-e','--prp_path', required=True)
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


def fakepair_corr_omic(mdata, omic_a, omic_b, barmap):
    x = mdata.mod[omic_a][barmap[omic_a.upper()]].X
    y = mdata.mod[omic_a][barmap[omic_b.upper()]].X
    df_cor = []
    obs_lst = mdata.obs_names
    var_lst = mdata[omic_a].var_names
    for i in range(x.shape[0]):
        obs_name = obs_lst[i]
        stat, pval = st.spearmanr(x[i, :].ravel(), y[i, :].ravel())
        df_cor.append(['obs', omic_a, obs_name, stat, pval])
    for j in range(x.shape[1]):
        var_name = var_lst[j]
        stat, pval = st.spearmanr(x[:, j].ravel(), y[:, j].ravel())
        df_cor.append(['var', omic_a, var_name, stat, pval])
    df_cor = pd.DataFrame(df_cor, columns=['type', 'omic', 'name', 'stat', 'pval'])
    return df_cor


df_k = []
df_cor = []
for pair in [('atac', 'rna'), ('rna', 'atac')]:
    omic_a, omic_b = pair
    df_k.append(calculate_k(mdata, omic_a, omic_b, barmap))
    df_cor.append(fakepair_corr_omic(mdata, omic_a, omic_b, barmap))
df_k = pd.concat(df_k)
df_cor = pd.concat(df_cor)

# Calculate confusion matrix
barmap.loc[:, 'ctype_rna'] = mdata[barmap['RNA']].obs['celltype'].values
barmap = barmap.rename(columns={'celltype': 'ctype_atac'})
cats_total = (
    barmap
    .groupby("ctype_atac")
    .size()
    .reset_index(name="n_total")
    .set_index("ctype_atac")
)

cats_df = (
    barmap
    .groupby(["ctype_atac", "ctype_rna"])
    .size()
    .reset_index(name="n")
    .set_index("ctype_atac")
    .join(cats_total)
    .reset_index()
)

cats_df = cats_df.assign(prop = lambda x: x.n / x.n_total)

# Write
df_k.to_csv(args.knn_path, index=False)
df_cor.to_csv(args.cor_path, index=False)
cats_df.to_csv(args.prp_path, index=False)
