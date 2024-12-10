import pandas as pd
import numpy as np
import mudata as mu
import decoupler as dc
import argparse
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils import f_beta_score


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--grn_path', required=True)
parser.add_argument('-p','--ptw_path', required=True)
parser.add_argument('-o','--out_path', required=True)
args = vars(parser.parse_args())

grn_path = args['grn_path']
ptw_path = args['ptw_path']
out_path = args['out_path']

grn_name = os.path.basename(grn_path).replace('.grn.csv', '')
data_path = os.path.join(os.path.dirname(os.path.dirname(grn_path)), 'mdata.h5mu')

grn = pd.read_csv(grn_path)

def get_sig_pws(grn, db, thr_pval):
    sig_pws = set()
    for tf in grn['source'].unique():
        df = grn[grn['source'] == tf].set_index('target')
        pws = dc.get_ora_df(
            df=df,
            net=db,
        )
        sig_pws.update(pws[pws['FDR p-value'] < thr_pval]['Term'])
    sig_pws = np.array(list(sig_pws))
    return sig_pws


def eval_metrics(y_pred, y):
    tp = np.intersect1d(y_pred, y).size
    if tp > 0.:
        fp = np.setdiff1d(y_pred, y).size
        fn = np.setdiff1d(y, y_pred).size
        prc = tp / (tp + fp)
        rcl = tp / (tp + fn)
        f1 = f_beta_score(prc, rcl)
    else:
        prc, rcl, f1 = 0., 0., 0.
    return prc, rcl, f1


def eval_grn(data, grn, db, thr_pval=0.01, thr_prop=0.2):
    hits = get_pw_hits(data, thr_pval, thr_prop)
    sig_pws = get_sig_pws(grn, db, thr_pval)
    prc, rcl, f1 = eval_metrics(y_pred=sig_pws, y=hits)
    return prc, rcl, f1


def get_pw_hits(data, thr_pval, thr_prop):
    pvals = data.obsm['ulm_pvals'].copy()
    pvals.loc[:, :] = dc.p_adjust_fdr(pvals.values.ravel()).reshape(pvals.shape)
    acts = data.obsm['ulm_estimate'].copy()
    hits = ((pvals < thr_pval) & (acts > 0)).sum(0).sort_values(ascending=False) / pvals.shape[0]
    hits = hits[hits > thr_prop].index.values.astype('U')
    return hits


if grn.shape[0] > 0:
    ptw = pd.read_csv(ptw_path)
    rna = mu.read(os.path.join(data_path, 'mod', 'rna'))
    # Infer pathway activities
    dc.run_ulm(
        mat=rna,
        net=ptw,
        weight=None,
        use_raw=False,
        verbose=True
    )
    prc, rcl, f01 = eval_grn(rna, grn, ptw, thr_pval=0.01, thr_prop=0.2)
    df = pd.DataFrame([[grn_name, prc, rcl, f01]], columns=['name', 'prc', 'rcl', 'f01'])
else:
    df = pd.DataFrame([[grn_name, np.nan, np.nan, np.nan]], columns=['name', 'prc', 'rcl', 'f01'])

# Write
df.to_csv(out_path, index=False)
