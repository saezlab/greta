from pyboolnet import file_exchange, trap_spaces
import scipy.stats as scs
import pandas as pd
import numpy as np
import decoupler as dc
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils import f_beta_score


def define_bool_rules(grn):
    targets = grn['target'].unique()
    rules = ''
    for target in targets:
        t_form = f'{target}, '
        msk = grn['target'] == target
        pos = ''
        neg = ''
        tgrn = grn.loc[msk]
        tgrn = tgrn.loc[tgrn['score'].abs().sort_values(ascending=False).index].head(10)
        for source, sign in zip(tgrn['source'], tgrn['score']):
            if sign >= 0:
                pos += f'{source} | '
            else:
                neg += f'!{source} & '
        pos, neg = pos[:-3], neg[:-3]
        if (pos != '') and (neg != ''):
            t_form += f'({pos}) & ({neg})'
        elif (pos != '') and (neg == ''):
            t_form += f'{pos}'
        elif (pos == '') and (neg != ''):
            t_form += f'{neg}'
        rules += t_form + '\n'
        print(t_form)
    return rules


def compute_fisher(hits, ct_set, tfs):
    a = hits & ct_set
    b = hits - ct_set
    c = ct_set - hits
    d = tfs - a - b - c
    a, b, c, d = len(a), len(b), len(c), len(d)
    p = dc.test1r(a, b, c, d)
    return p


def find_hits(sss, ct_sets, tfs, thr_pval):
    pvals = np.zeros((len(sss), ct_sets.shape[0]))
    for i, pss in enumerate(sss):
        hits = set()
        for k in pss:
            if pss[k]:
                hits.add(k)
        for j, ct_set in enumerate(ct_sets):
            pvals[i, j] = compute_fisher(hits, ct_set, tfs)
    pvals = scs.false_discovery_control(pvals, axis=1)
    df = pd.DataFrame(pvals < thr_pval, columns=ct_sets.index)
    return df


def get_prc_rcl(hits):
    n_ct_per_ss = hits.sum(1)
    tp = (n_ct_per_ss == 1).sum()
    if tp > 0:
        fp = (n_ct_per_ss != 1).sum()
        prc = tp / (tp + fp)
    else:
        prc = 0
    
    n_cts = hits.sum(0)
    tp = (n_cts > 0).sum()
    if tp > 0:
        tn = (n_cts == 0).sum()
        rcl = tp / (tp + tn)
    else:
        rcl = 0
    return prc, rcl


def compute_score(grn, ct_df, thr_pval=0.01):
    # Find ct sets
    ct_sets = ct_df.groupby('celltype')['tf'].apply(lambda x: set(x))
    # Filter for tfs
    tfs = set(ct_df['tf'])
    sgrn = grn[(grn['source'].isin(tfs)) & (grn['target'].isin(tfs))]
    if sgrn.shape[0] >= 5:
        # Simulate steady states
        rules = define_bool_rules(sgrn)
        print(rules)
        print('Generating primes ...')
        primes = file_exchange.bnet2primes(rules)
        print('Primes generated')
        print('Computing steady states ...')
        sss = trap_spaces.compute_steady_states(primes, max_output=int(100_000))
        hits = find_hits(sss, ct_sets, tfs, thr_pval)
        print('Done')
        prc, rcl = get_prc_rcl(hits)
        f01 = f_beta_score(prc, rcl)
        return prc, rcl, f01
    else:
        return np.nan, np.nan, np.nan


path_grn = sys.argv[1]
path_tfs = sys.argv[2]
thr_pval = float(sys.argv[3])
path_out = sys.argv[4]

# Load data
grn_name = os.path.basename(path_grn).replace('.grn.csv', '')
grn = pd.read_csv(path_grn).drop_duplicates(['source', 'target'])
ct_df = pd.read_csv(path_tfs)

# Compute score
prc, rcl, f01 = compute_score(grn, ct_df)

# Transform to df
df = pd.DataFrame([[grn_name, prc, rcl, f01]], columns=['name', 'prc', 'rcl', 'f01'])

# Write
df.to_csv(path_out, index=False)
