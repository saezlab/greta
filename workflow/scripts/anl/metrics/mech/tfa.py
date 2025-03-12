import pandas as pd
import numpy as np
import decoupler as dc
import mudata as mu
import sys
import os
import re
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils import load_cats, f_beta_score
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

# Extract names and path
grn_name = os.path.basename(grn_path).replace('.grn.csv', '')
data_path = os.path.join(os.path.dirname(os.path.dirname(grn_path)), 'mdata.h5mu')
dataset = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(data_path))))
case = os.path.basename(os.path.dirname(data_path))
rsc_name = os.path.basename(bnc_path)

# Read GRN
grn = pd.read_csv(grn_path)
grn = grn.drop_duplicates(['source', 'target'], keep='first')

if grn.shape[0] > 0:
    # Read dataset
    rna = mu.read(os.path.join(data_path, 'mod', 'rna'))

    # Read benchmark data
    mat = pd.read_csv(os.path.join(bnc_path, 'diff.csv'), index_col=0)
    obs = pd.read_csv(os.path.join(bnc_path, 'meta.csv'), index_col=0)
    
    # Subset bench data to dataset
    cats = load_cats(dataset, case)
    cats = [re.escape(c) for c in cats[rsc_name]]
    msk = obs['Tissue.Type'].isin(cats) & obs['TF'].isin(rna.var_names) & (obs['logFC'] < -0.5)
    obs = obs.loc[msk, :]
    mat = mat.loc[msk, :]

    # Compute TF activities
    acts = []
    pvals = []
    for dataset in obs.index:
        tf = obs.loc[dataset, 'TF']
        tf_mat = mat.loc[[dataset], :]
        tf_grn = grn[grn['source'] == tf]
        try:
            act, pval = dc.run_ulm(
                mat=tf_mat,
                net=tf_grn,
                weight='score',
                min_n=3,
            )
            act, pval = act.values[0, 0], pval.values[0, 0]
            acts.append(act)
            pvals.append(pval)
        except:
            pass

    # Compute recall
    acts = np.array(acts)
    pvals = np.array(pvals)
    padj = dc.p_adjust_fdr(pvals)
    tp = np.sum((acts < 0) & (padj < 0.05))
    if tp > 0:
        prc = tp / acts.size
        rcl = tp / obs.shape[0]
        f01 = f_beta_score(prc, rcl)
    else:
        prc, rcl, f01 = 0., 0., 0.

    df = pd.DataFrame([[grn_name, prc, rcl, f01]], columns=['name', 'prc', 'rcl', 'f01'])
else:
    df = pd.DataFrame([[grn_name, np.nan, np.nan, np.nan]], columns=['name', 'prc', 'rcl', 'f01'])

# Write
df.to_csv(out_path, index=False)
