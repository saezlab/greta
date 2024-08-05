import pandas as pd
import numpy as np
import decoupler as dc
import mudata as mu
import celloracle as co
import scipy
import os
from tqdm import tqdm
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--grn_path', required=True)
parser.add_argument('-b','--bnc_path', required=True)
parser.add_argument('-c','--cats', required=True, nargs='+')
parser.add_argument('-o','--out_path', required=True)
args = vars(parser.parse_args())

grn_path = args['grn_path']
bnc_path = args['bnc_path']
cats = args['cats']
out_path = args['out_path']

# Extract names and path
grn_name = os.path.basename(grn_path).replace('.grn.csv', '')
data_path = os.path.join(os.path.dirname(os.path.dirname(grn_path)), 'mdata.h5mu')

# Read GRN
grn = pd.read_csv(grn_path)
grn = grn.drop_duplicates(['source', 'target'], keep='first')

def init_celloracle(adata, grn):
    oracle = co.Oracle()
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
    col_dict = co.trajectory.oracle_utility._get_clustercolor_from_anndata(adata=oracle.adata,
                                              cluster_name='cluster',
                                              return_as="dict")
    oracle.colorandum = np.array([col_dict[i] for i in oracle.adata.obs['cluster']])
    # Add grn
    tf_dict = grn.groupby(['target'])['source'].apply(lambda x: sorted(list(x))).to_dict()   # Refit GRN
    oracle.addTFinfo_dictionary(tf_dict)
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
    
    # Read benchmark data
    mat = pd.read_csv(os.path.join(bnc_path, 'log2fcs.csv'), index_col=0)
    obs = pd.read_csv(os.path.join(bnc_path, 'meta.csv'), index_col=0)
    
    # Subset bench data to dataset
    msk = obs['Tissue.Type'].isin(cats) & obs['TF'].isin(rna.var_names) & (obs['logFC'] < -0.5)
    obs = obs.loc[msk, :]
    mat = mat.loc[msk, :]

    oracle = init_celloracle(rna, grn)
    coefs = []
    pvals = []
    for dataset in tqdm(obs.index):
        # Extract
        tf = obs.loc[dataset, 'TF']
        tf_mat = mat.loc[[dataset], :]
    
        try:
            # Run the simulation for the current TF
            delta = simulate_delta(oracle, [tf], n_steps=3)
            
            # Compute correlation
            x = delta.mean(0)
            y = tf_mat.T[dataset]
            inter = np.intersect1d(x.index, y.index)
            x, y = x.loc[inter].values, y.loc[inter]
            r, p = scipy.stats.spearmanr(x, y)
            coefs.append(r)
            pvals.append(p)
        except:
            pass
    
    # Compute recall
    coefs = np.array(coefs)
    pvals = np.array(pvals)
    padj = dc.p_adjust_fdr(pvals)
    tp = np.sum((coefs > 0.05) & (padj < 0.05))
    if tp > 0:
        prc = tp / coefs.size
        rcl = tp / obs.shape[0]
        beta = 0.1
        f01 = ((1 - beta**2) * prc * rcl) / ((prc * beta**2) + rcl)
    else:
        prc, rcl, f01 = 0., 0., 0.
    df = pd.DataFrame([[grn_name, prc, rcl, f01]], columns=['name', 'prc', 'rcl', 'f01'])
else:
    df = pd.DataFrame([[grn_name, np.nan, np.nan, np.nan]], columns=['name', 'prc', 'rcl', 'f01'])

# Write
df.to_csv(out_path, index=False)
