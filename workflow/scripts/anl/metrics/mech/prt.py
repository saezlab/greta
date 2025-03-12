import pandas as pd
import numpy as np
import mudata as mu
import anndata as ad
import sys
import os
import re
import appdirs
import tempfile
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--grn_path', required=True)
parser.add_argument('-b','--bnc_path', required=True)
parser.add_argument('-o','--out_path', required=True)
args = vars(parser.parse_args())

grn_path = args['grn_path']
bnc_path = args['bnc_path']
out_path = args['out_path']

grn_name = os.path.basename(grn_path).replace('.grn.csv', '')
def user_cache_dir(*args, **kwargs):
    tmp_dir = tempfile.gettempdir()
    return os.path.join(tmp_dir, grn_name)
sys.modules["appdirs"].user_cache_dir = user_cache_dir

from celloracle import Oracle
import celloracle.trajectory.oracle_utility as co
import scipy
from tqdm import tqdm
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils import load_cats, f_beta_score

# Extract names and path
data_path = os.path.join(os.path.dirname(os.path.dirname(grn_path)), 'mdata.h5mu')
dataset = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(data_path))))
case = os.path.basename(os.path.dirname(data_path))
rsc_name = os.path.basename(bnc_path)

# Read GRN
grn = pd.read_csv(grn_path)
grn = grn.drop_duplicates(['source', 'target'], keep='first')

def init_celloracle(adata, grn, fit_grn):
    oracle = Oracle()
    oracle.adata = adata
    oracle.adata.obsm['X_umap'] = np.zeros((adata.shape[0], 2))
    oracle.adata.layers['imputed_count'] = oracle.adata.X
    oracle.adata.obs['cluster'] = 'cluster'
    oracle.cluster_column_name = None
    oracle.embedding_name = 'X_umap'
    oracle.pcs = np.zeros((oracle.adata.shape[0], 2))
    oracle.knn = True
    oracle.k_knn_imputation = True
    oracle.high_var_genes = list(oracle.adata.var_names)
    oracle.adata.obs['cluster'] = oracle.adata.obs['cluster'].astype('category')
    oracle.adata.uns['cluster_colors'] = ['#1f77b4']
    col_dict = co._get_clustercolor_from_anndata(adata=oracle.adata,
                                              cluster_name='cluster',
                                              return_as="dict")
    oracle.colorandum = np.array([col_dict[i] for i in oracle.adata.obs['cluster']])
    tf_dict = grn.groupby(['target'])['source'].apply(lambda x: sorted(list(x))).to_dict()   # Refit GRN
    oracle.addTFinfo_dictionary(tf_dict)
    # Add grn
    if fit_grn:
        oracle.fit_GRN_for_simulation(alpha=0, GRN_unit="whole")
    return oracle


def simulate_delta(oracle, tfs, n_steps=3):
    perturb_condition = dict()

    # Knockdown the specified transcription factors
    for g in tfs:
        perturb_condition[g] = 0  # Use zero for knockdown

    # Perform the simulation
    oracle.simulate_shift(perturb_condition=perturb_condition, n_propagation=n_steps)

    # Retrieve and process the delta values
    delta = pd.DataFrame(oracle.adata.layers['delta_X'], index=oracle.adata.obs_names, columns=oracle.adata.var_names)

    # Remove the columns that are 0s.
    delta = delta.loc[:, delta.abs().sum(0) != 0]

    return delta


if grn.shape[0] > 0:
    # Read dataset
    rna = mu.read(os.path.join(data_path, 'mod', 'rna'))

    # Subset data to grn
    genes = set(grn['source']) | set(grn['target'])
    rna = rna[:, rna.var_names.isin(genes)]

    # Init object
    oracle = init_celloracle(rna, grn, fit_grn=True)
    coef_mat = oracle.coef_matrix
    oracle = ad.AnnData(oracle.adata.to_df().mean(0).to_frame().T)
    oracle = init_celloracle(oracle, grn, fit_grn=False)
    oracle.coef_matrix = coef_mat
    tf_n_trgs = (coef_mat != 0).sum(0)
    tf_n_trgs = set(tf_n_trgs[tf_n_trgs >= 3].index)
    oracle.all_regulatory_genes_in_TFdict = [t for t in oracle.all_regulatory_genes_in_TFdict if t in tf_n_trgs]

    # Read benchmark data
    mat = pd.read_csv(os.path.join(bnc_path, 'diff.csv'), index_col=0)
    obs = pd.read_csv(os.path.join(bnc_path, 'meta.csv'), index_col=0)

    # Subset bench data to dataset
    cats = load_cats(dataset, case)
    cats = [re.escape(c) for c in cats[rsc_name]]
    msk = obs['Tissue.Type'].isin(cats) & obs['TF'].isin(rna.var_names) & (obs['logFC'] < -0.5)
    obs = obs.loc[msk, :]
    mat = mat.loc[msk, :]

    # Subset by overlap with rna
    genes = list(genes & set(mat.columns))
    mat = mat.loc[:, genes].copy()
    
    coefs = []
    pvals = []
    for dataset in tqdm(obs.index):
        # Extract
        tf = obs.loc[dataset, 'TF']
        tf_mat = mat.loc[[dataset], :]
        tf_mat = tf_mat[tf_mat != 0].dropna(axis=1)
        if tf in oracle.all_regulatory_genes_in_TFdict:
            # Run the simulation for the current TF
            x = simulate_delta(oracle, [tf], n_steps=3)
            
            # Intersect
            y = tf_mat
            inter = np.intersect1d(x.columns, y.columns)
            x, y = x.loc[:, inter].values[0], y.loc[:, inter].values[0]

            # Compute correlation
            if x.size >= 10:
                r, p = scipy.stats.spearmanr(x, y)
            else:
                r, p = 0., 1.
            coefs.append(r)
            pvals.append(p)
    
    # Compute recall
    coefs = np.array(coefs)
    pvals = np.array(pvals)
    padj = scipy.stats.false_discovery_control(pvals, method='bh')
    tp = np.sum((coefs > 0.05) & (padj < 0.05))
    if tp > 0:
        prc = tp / coefs.size
        rcl = tp / obs.shape[0]
        f01 = f_beta_score(prc, rcl)
    else:
        prc, rcl, f01 = 0., 0., 0.
    df = pd.DataFrame([[grn_name, prc, rcl, f01]], columns=['name', 'prc', 'rcl', 'f01'])
else:
    df = pd.DataFrame([[grn_name, np.nan, np.nan, np.nan]], columns=['name', 'prc', 'rcl', 'f01'])

# Write
df.to_csv(out_path, index=False)
