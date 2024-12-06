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


# Read
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
    ctype_dict = (
        mdata.obs
        .reset_index(names='barcode')
        .groupby('celltype')['barcode']
        .apply(lambda x: list(x))
    )
    ctypes = mdata.obs['celltype'].values
    # Find closest k
    rng = np.random.default_rng(seed=42)
    for barcode_a, barcode_b in zip(barmap[omic_a.upper()], barmap[omic_b.upper()]):
        i = bar_lst.index(barcode_a)
        ctyp = mdata.obs.loc[barcode_a, 'celltype']
        knn = indices[i]
        sorted_barcodes = list(omic.obs_names.values[knn])
        k = sorted_barcodes.index(barcode_b) + 1
        df.append(['predicted', ctyp, omic_a, barcode_a, k])
        barcode_b = rng.choice(ctype_dict[ctypes[i]], size=1)[0]
        ctyp = mdata.obs.loc[barcode_b, 'celltype']
        k = sorted_barcodes.index(barcode_b) + 1
        df.append(['random', ctyp, omic_a, barcode_a, k])
    df = pd.DataFrame(df, columns=['type', 'ctype', 'anchor', 'barcode', 'k'])
    return df


def fakepair_corr_omic(mdata, omic_a, omic_b, barmap):
    x = mdata.mod[omic_a][barmap[omic_a.upper()]].X
    y = mdata.mod[omic_a][barmap[omic_b.upper()]].X
    ctype_dict = (
        mdata[barmap[omic_a.upper()]].obs
        .reset_index(names='barcode')
        .reset_index()
        .groupby('celltype')['index']
        .apply(lambda x: list(x))
    )
    ctypes = mdata.obs['celltype'].values
    df_cor = []
    obs_lst = mdata.obs_names
    var_lst = mdata[omic_a].var_names
    rng = np.random.default_rng(seed=42)
    for i in range(x.shape[0]):
        obs_name = obs_lst[i]
        ctyp = ctypes[i]
        stat, pval = st.spearmanr(x[i, :].ravel(), y[i, :].ravel())
        df_cor.append(['predicted', ctyp, omic_a, obs_name, stat, pval])
        r = rng.choice(ctype_dict[ctyp], size=1)[0]
        stat, pval = st.spearmanr(x[i, :].ravel(), x[r, :].ravel())
        df_cor.append(['random', ctyp, omic_a, obs_name, stat, pval])
    df_cor = pd.DataFrame(df_cor, columns=['type', 'ctype', 'omic', 'name', 'stat', 'pval'])
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
