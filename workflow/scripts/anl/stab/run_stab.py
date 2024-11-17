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
parser.add_argument('-o','--out_path', required=True)
args = vars(parser.parse_args())

inp_path = args['inp_path']
out_path = args['out_path']

def time_to_hours(time_str):
    h, m, s = map(int, time_str.split(':'))
    return h + (m / 60) + (s / 60 / 60)


def memory_to_gb(memory_str):
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
dataset = os.path.basename(out_path).replace('.csv', '')
res = []
for _, row in list(df.iterrows()):
    s, time, mem = row[0], row[1], row[2]
    mth = re.search(r'mdl_o_(.*?)_dat', s).group(1)
    ds = re.search(r'=(.*?).case=', s).group(1)
    case = re.search(r'case=([\w_]+)', s).group(1)
    
    if (ds == dataset) or (case != 'all'):
        if case.startswith('16384'):
            if case.startswith('16384_16384_'):
                cat = 'full'
                ncells, nfeats, seed = case.split('_')
                n = 16384
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
        ref = pd.read_csv('dts/{dataset}/cases/16384_16384_0/runs/o_{mth}.o_{mth}.o_{mth}.o_{mth}.grn.csv'.
                          format(dataset=ds, mth=mth))
        net = pd.read_csv('dts/{dataset}/cases/{ncells}_{nfeats}_{seed}/runs/o_{mth}.o_{mth}.o_{mth}.o_{mth}.grn.csv'.
                          format(dataset=ds, ncells=ncells, nfeats=nfeats, seed=seed, mth=mth))
        tmp = pd.DataFrame(index=[0])
        tmp['mth'] = mth
        tmp['cat'] = cat
        tmp['n'] = int(n)
        tmp['seed'] = int(seed)
        tmp['h'] = time_to_hours(time)
        tmp['gb'] = memory_to_gb(mem)
        tmp['ocoeff'] = ocoeff(ref, net)
        tmp[['n_sources', 'n_edges', 'n_targets', 'r_size']] = get_grn_stats(net)
        res.append(tmp)
res = pd.concat(res)

# Drop potential duplicates from sacct
res = res.drop_duplicates(['mth', 'cat', 'n', 'seed'], keep='last')

# Sort
res = res.sort_values(['mth', 'cat', 'n', 'seed']).reset_index(drop=True)

# Write
res.to_csv(out_path, index=False)
