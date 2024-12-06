import pandas as pd
import numpy as np
import argparse
import os
import re
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils import ocoeff


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--inp_path', required=True)
parser.add_argument('-r','--res_path', required=True)
parser.add_argument('-a','--auc_path', required=True)
args = vars(parser.parse_args())

inp_path = args['inp_path']
res_path = args['res_path']
auc_path = args['auc_path']

def time_to_hours(time_str):
    h, m, s = map(int, time_str.split(':'))
    return h + (m / 60) + (s / 60 / 60)


def memory_to_gb(memory_str):
    if memory_str != '0':
        unit = memory_str[-1]
        value = int(np.ceil(float(memory_str[:-1])))
        if unit == 'K':
            return value / 1048576.0
        elif unit == 'M':
            return value / 1024.0
        elif unit == 'G':
            return float(value)
        else:
            raise ValueError("Unsupported memory unit. Please use K, M, or G.")
    else:
        return 0.


def get_grn_stats(df):
    n_s = df['source'].unique().size
    n_e = df.shape[0]
    n_t = df['target'].unique().size
    n_r = df.groupby(['source']).count()['target'].mean()
    if np.isnan(n_r):
        n_r = 0.
    return n_s, n_e, n_t, n_r


# Read
df = pd.read_csv(inp_path, sep=' ', header=None)
dataset = os.path.basename(res_path).replace('.csv', '')
res = []
seeds = np.array(['0', '1', '2'])
for _, row in list(df.iterrows()):
    s, time, mem = row[0], row[1], row[2]
    mth = re.search(r'mdl_(.*?)_dat', s).group(1)
    ds = re.search(r'=(.*?).case=', s).group(1)
    case = re.search(r'case=([\w_]+)', s).group(1)
    
    if (ds == dataset) or (case != 'all'):
        if case.startswith('16384'):
            if case.startswith('16384_16384_'):
                cat = 'full'
                ncells, nfeats, seed = case.split('_')
                n = 16384
                net = pd.read_csv('dts/{dataset}/cases/{ncells}_{nfeats}_{seed}/runs/{mth}.{mth}.{mth}.{mth}.grn.csv'.
                          format(dataset=ds, ncells=ncells, nfeats=nfeats, seed=seed, mth=mth))
                for s in [s for s in seeds if s != seed]:
                    ref = pd.read_csv('dts/{dataset}/cases/{ncells}_{nfeats}_{seed}/runs/{mth}.{mth}.{mth}.{mth}.grn.csv'.
                          format(dataset=ds, ncells=ncells, nfeats=nfeats, seed=s, mth=mth))
                    tmp = pd.DataFrame(index=[0])
                    tmp['mth'] = mth.replace('o_', '')
                    tmp['cat'] = cat
                    tmp['n'] = int(n)
                    tmp['seed'] = int(seed)
                    tmp['other_seed'] = int(s)
                    tmp['h'] = time_to_hours(time)
                    tmp['gb'] = memory_to_gb(mem)
                    tmp['s_ocoeff'] = ocoeff(ref, net, on=['source'])
                    tmp['e_ocoeff'] = ocoeff(ref, net, on=['source', 'target'])
                    tmp['t_ocoeff'] = ocoeff(ref, net, on=['target'])
                    tmp[['n_sources', 'n_edges', 'n_targets', 'r_size']] = get_grn_stats(net)
                    res.append(tmp)
                continue
            else:
                cat = 'fixed_ncells'
                _, n, seed = case.split('_')
                ncells, nfeats = 16384, n
        elif '_16384_' in case:
            cat = 'fixed_nfeats'
            n, _, seed = case.split('_')
            ncells, nfeats = n, 16384
        else:
            continue
        ref = pd.read_csv('dts/{dataset}/cases/16384_16384_0/runs/{mth}.{mth}.{mth}.{mth}.grn.csv'.
                          format(dataset=ds, mth=mth))
        net = pd.read_csv('dts/{dataset}/cases/{ncells}_{nfeats}_{seed}/runs/{mth}.{mth}.{mth}.{mth}.grn.csv'.
                          format(dataset=ds, ncells=ncells, nfeats=nfeats, seed=seed, mth=mth))
        tmp = pd.DataFrame(index=[0])
        tmp['mth'] = mth.replace('o_', '')
        tmp['cat'] = cat
        tmp['n'] = int(n)
        tmp['seed'] = int(seed)
        tmp['other_seed'] = np.nan
        tmp['h'] = time_to_hours(time)
        tmp['gb'] = memory_to_gb(mem)
        tmp['s_ocoeff'] = ocoeff(ref, net, on=['source'])
        tmp['e_ocoeff'] = ocoeff(ref, net, on=['source', 'target'])
        tmp['t_ocoeff'] = ocoeff(ref, net, on=['target'])
        tmp[['n_sources', 'n_edges', 'n_targets', 'r_size']] = get_grn_stats(net)
        res.append(tmp)
res = pd.concat(res)

# Drop potential duplicates from sacct
res = res.drop_duplicates(['mth', 'cat', 'n', 'seed', 'other_seed'], keep='last')

# Sort
res = res.sort_values(['mth', 'cat', 'n', 'seed']).reset_index(drop=True)

# Compute stab score (auc)
mdf = res.groupby(['mth', 'cat', 'n']).mean().reset_index()
aucs = []
types = ['s_ocoeff', 'e_ocoeff', 't_ocoeff']
for mth in mdf['mth'].sort_values().unique():
    for cat in ['fixed_nfeats', 'fixed_ncells']:
        tmp = mdf[(mdf['mth'] == mth) & (mdf['cat'].isin([cat, 'full']))]
        for typ in types:
            y = tmp[typ].values
            x = np.arange(y.size) / (y.size - 1)
            auc = np.trapezoid(y, x)
            aucs.append([typ, mth, cat, auc])
acus = pd.DataFrame(aucs, columns=['type', 'mth', 'cat', 'auc'])

# Write
res.to_csv(res_path, index=False)
acus.to_csv(auc_path, index=False)
