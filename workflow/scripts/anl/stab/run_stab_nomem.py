import pandas as pd
import numpy as np
import argparse
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils import ocoeff


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--inp_paths', nargs='+', required=True)
parser.add_argument('-n', '--n_seeds', required=True)
parser.add_argument('-r', '--res_path', required=True)
parser.add_argument('-a', '--auc_path', required=True)
args = vars(parser.parse_args())

inp_paths = args['inp_paths']
n_seeds = int(args['n_seeds'])
res_path = args['res_path']
auc_path = args['auc_path']


_cache = {}


def read_grn(path):
    if path not in _cache:
        _cache[path] = pd.read_csv(path)
    return _cache[path]


# Parse all paths
# e.g. dts/hg38/pitupair/cases/16384_1024_3/runs/o_celloracle.o_celloracle.o_celloracle.o_celloracle.grn.csv
parsed = []
for path in inp_paths:
    parts = path.split('/')
    org = parts[1]
    ds = parts[2]
    case = parts[4]
    mth = parts[6].split('.')[0]
    ncells, nfeats, seed = case.split('_')
    parsed.append({'path': path, 'org': org, 'ds': ds, 'mth': mth,
                   'ncells': ncells, 'nfeats': nfeats, 'seed': seed})

# Group by (org, ds, mth) for reference lookups
groups = defaultdict(list)
for p in parsed:
    groups[(p['org'], p['ds'], p['mth'])].append(p)

res = []

for (org, ds, mth), items in groups.items():
    # Index full (16384_16384_*) paths by seed for reference lookups
    full_items = {
        item['seed']: item
        for item in items
        if item['ncells'] == '16384' and item['nfeats'] == '16384'
    }

    for item in items:
        ncells, nfeats, seed = item['ncells'], item['nfeats'], item['seed']

        if ncells == '16384' and nfeats == '16384':
            cat = 'full'
            n = 16384
            net = read_grn(item['path'])
            for other_seed, other_item in full_items.items():
                if other_seed == seed:
                    continue
                ref = read_grn(other_item['path'])
                tmp = pd.DataFrame(index=[0])
                tmp['mth'] = mth.replace('o_', '')
                tmp['cat'] = cat
                tmp['n'] = int(n)
                tmp['seed'] = int(seed)
                tmp['other_seed'] = int(other_seed)
                tmp['s_ocoeff'] = ocoeff(ref, net, on=['source'])
                tmp['c_ocoeff'] = ocoeff(ref, net, on=['cre'])
                tmp['t_ocoeff'] = ocoeff(ref, net, on=['target'])
                tmp['e_ocoeff'] = ocoeff(ref, net, on=['source', 'target'])
                res.append(tmp)

        elif ncells == '16384':
            cat = 'fixed_ncells'
            n = nfeats
            if '0' not in full_items:
                continue
            ref = read_grn(full_items['0']['path'])
            net = read_grn(item['path'])
            tmp = pd.DataFrame(index=[0])
            tmp['mth'] = mth.replace('o_', '')
            tmp['cat'] = cat
            tmp['n'] = int(n)
            tmp['seed'] = int(seed)
            tmp['other_seed'] = np.nan
            tmp['s_ocoeff'] = ocoeff(ref, net, on=['source'])
            tmp['c_ocoeff'] = ocoeff(ref, net, on=['cre'])
            tmp['t_ocoeff'] = ocoeff(ref, net, on=['target'])
            tmp['e_ocoeff'] = ocoeff(ref, net, on=['source', 'target'])
            res.append(tmp)

        elif nfeats == '16384':
            cat = 'fixed_nfeats'
            n = ncells
            if '0' not in full_items:
                continue
            ref = read_grn(full_items['0']['path'])
            net = read_grn(item['path'])
            tmp = pd.DataFrame(index=[0])
            tmp['mth'] = mth.replace('o_', '')
            tmp['cat'] = cat
            tmp['n'] = int(n)
            tmp['seed'] = int(seed)
            tmp['other_seed'] = np.nan
            tmp['s_ocoeff'] = ocoeff(ref, net, on=['source'])
            tmp['c_ocoeff'] = ocoeff(ref, net, on=['cre'])
            tmp['t_ocoeff'] = ocoeff(ref, net, on=['target'])
            tmp['e_ocoeff'] = ocoeff(ref, net, on=['source', 'target'])
            res.append(tmp)

        else:
            continue

res = pd.concat(res)

# Drop potential duplicates and sort
res = res.drop_duplicates(['mth', 'cat', 'n', 'seed', 'other_seed'], keep='last')
res = res.sort_values(['mth', 'cat', 'n', 'seed']).reset_index(drop=True)

# Compute stab score (auc)
mdf = res.groupby(['mth', 'cat', 'n']).mean().reset_index()
aucs = []
types = ['s_ocoeff', 'c_ocoeff', 't_ocoeff', 'e_ocoeff']
for mth_name in mdf['mth'].sort_values().unique():
    for cat in ['fixed_nfeats', 'fixed_ncells']:
        tmp = mdf[(mdf['mth'] == mth_name) & (mdf['cat'].isin([cat, 'full']))]
        for typ in types:
            y = tmp[typ].values
            x = np.arange(y.size) / (y.size - 1)
            auc = np.trapezoid(y, x)
            aucs.append([typ, mth_name, cat, auc])
acus = pd.DataFrame(aucs, columns=['type', 'mth', 'cat', 'auc'])

# Write
res.to_csv(res_path, index=False)
acus.to_csv(auc_path, index=False)
