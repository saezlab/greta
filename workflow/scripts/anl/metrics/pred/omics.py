import pandas as pd
import numpy as np
import mudata as mu
from sklearn.model_selection import train_test_split
from xgboost import XGBRegressor
from tqdm import tqdm
import scipy
import argparse
import sys
import pyranges as pr
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils import f_beta_score


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-a','--grn_path', required=True)
parser.add_argument('-b','--col_source', required=True)
parser.add_argument('-c','--col_target', required=True)
parser.add_argument('-d','--mod_source', required=True)
parser.add_argument('-e','--mod_target', required=True)
parser.add_argument('-f','--out_path', required=True)
args = vars(parser.parse_args())

grn_path = args['grn_path']
col_source = args['col_source']
col_target = args['col_target']
mod_source = args['mod_source']
mod_target = args['mod_target']
out_path = args['out_path']

grn_name = os.path.basename(grn_path).replace('.grn.csv', '')
data_path = os.path.join(os.path.dirname(os.path.dirname(grn_path)), 'mdata.h5mu')

grn = pd.read_csv(grn_path)


def map_regions(a, b):
    def to_df(regions):
        data = [r.split('-') for r in regions]
        return pd.DataFrame(data, columns=['Chromosome','Start','End']).assign(Start=lambda x: x.Start.astype(int), End=lambda x: x.End.astype(int), region=regions)
    pr_a = pr.PyRanges(to_df(a))
    pr_b = pr.PyRanges(to_df(b))
    joined = pr_b.join(pr_a, suffix='_a')
    return joined.df[['region','region_a']].rename(columns={'region':'cre','region_a':'new_cre'})


def test_predictability(mdata, train, test, grn, col_source='source', col_target='target', mod_source='rna', mod_target='rna', ntop=5):
    def remove_zeros(X, y):
        msk = y != 0.
        y = y[msk]
        X = X[msk, :]
        return X, y
    net = grn.drop_duplicates([col_source, col_target]).copy()
    net['abs_score'] = abs(net['score'])
    net = net.sort_values('abs_score', ascending=False)
    net = net.groupby(col_target)[col_source].apply(lambda x: list(x) if ntop is None else list(x)[:ntop])
    cor = []
    for target in tqdm(net.index):
        sources = net[target]
        sources = [s for s in sources if s != target]
        if len(sources) > 0:
            train_X = mdata.mod[mod_source][train, sources].X
            train_y = mdata.mod[mod_target][train, target].X.ravel()
            test_X = mdata.mod[mod_source][test, sources].X
            test_y = mdata.mod[mod_target][test, target].X.ravel()
            train_X, train_y = remove_zeros(train_X, train_y)
            test_X, test_y = remove_zeros(test_X, test_y)
            if test_y.size >= 10:
                reg = XGBRegressor(random_state=0, n_jobs=1).fit(train_X, train_y)
                pred_y = reg.predict(test_X)
                if np.any(pred_y != pred_y[0]):
                    s, p = scipy.stats.spearmanr(pred_y, test_y)  # Spearman to control for outliers
                    cor.append([target, pred_y.size, len(sources), s, p])
    cor = pd.DataFrame(cor, columns=['target', 'n_obs', 'n_vars', 'coef', 'pval'])
    if cor.shape[0] > 0:
        cor['padj'] = scipy.stats.false_discovery_control(cor['pval'], method='bh')
    else:
        cor['padj'] = pd.Series(dtype=float)
    return cor

if grn.shape[0] > 0:
    mdata = mu.read_h5mu(data_path)
    # Filter GRN
    genes = mdata.mod['rna'].var_names
    peaks = mdata.mod['atac'].var_names
    cres = grn['cre'].unique()
    mapped_regions = map_regions(list(peaks), list(cres))
    grn = pd.merge(grn, mapped_regions, on='cre', how='inner').drop(columns=['cre']).rename(columns={'new_cre': 'cre'})
    grn = grn[(grn['target'].isin(genes)) & (grn['source'].isin(genes))]
    # Split and test
    train, test = train_test_split(mdata.obs_names, test_size=0.33, random_state=42, stratify=mdata.obs['celltype'])
    cor = test_predictability(mdata=mdata, train=train, test=test, grn=grn, col_source=col_source, col_target=col_target, mod_source=mod_source, mod_target=mod_target)
    sig_cor = cor[(cor['padj'] < 0.05) & (cor['coef'] > 0.05)]
    n_hits = sig_cor.shape[0]
    if n_hits > 0:
        universe_size = mdata.mod[mod_target].var_names.size
        rcl = n_hits / universe_size
        prc = n_hits / cor.shape[0]
        f01 = f_beta_score(prc, rcl)
    else:
        prc, rcl, f01 = 0., 0., 0.
    df = pd.DataFrame([[grn_name, prc, rcl, f01]], columns=['name', 'prc', 'rcl', 'f01'])
else:
    df = pd.DataFrame([[grn_name, np.nan, np.nan, np.nan]], columns=['name', 'prc', 'rcl', 'f01'])

# Write
df.to_csv(out_path, index=False)
