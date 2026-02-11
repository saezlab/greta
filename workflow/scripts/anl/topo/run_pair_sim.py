import concurrent.futures
import pandas as pd
import numpy as np
import os
import glob
from tqdm import tqdm
from functools import partial
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils import (
    ocoeff,
    get_grn_name,
    get_grn_stats,
    _map_regions,
)
import argparse
import pyranges as pr


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-t','--stat_path', required=True)
parser.add_argument('-s','--sim_path', required=True)
args = vars(parser.parse_args())

stat_path = args['stat_path']
sim_path = args['sim_path']

org, dat, case = os.path.basename(stat_path).split('.')[:3]
paths = glob.glob(os.path.join('dts', org, dat, 'cases', case, 'runs', '*.grn.csv'))

print('Reading and computing grns stats...')
names = []
dfs = []
stats = []
tfs = []
cres = []
genes = []
edges = []

for path in tqdm(paths):
    name = get_grn_name(path)
    names.append(name)
    df = pd.read_csv(path)
    stat = get_grn_stats(df)
    stats.append([name] + list(stat))
    tfs.append(set(df['source']))
    cres.append(set(df['cre']))
    genes.append(set(df['target']))
    tdf = df.drop_duplicates(['source', 'target'], keep='first')
    edges.append(set(tdf['source'] + '|' + tdf['target']))


# Store as df
cols = ['name', 'n_tfs', 'n_cres', 'n_targets', 'n_edges', 'odegree', 'betweenc', 'eigv']
stats = pd.DataFrame(stats, columns=cols)

print('Computing pairwise overlap coefficients...')


#def set_ocoef(a, b):
#    min_s = min(len(a), len(b))
#    if min_s == 0:
#        return np.nan
#    else:
#        inter = len(a & b)
#        return inter / min_s

def set_ocoef(a, b, use_overlap: bool = False) -> float:
    min_s = min(len(a), len(b))
    a_size = len(a)
    b_size = len(b)
    if min_s > 0:
        i_size = len(a & b)
        coeff = i_size / min_s
        if coeff == 0 and use_overlap:
            mapping = _map_regions(list(a), list(b))
            if mapping.empty:
                i_size = 0
            else:
                # Count overlapping regions from the smaller set
                if a_size <= b_size:
                    i_size = mapping["region_a"].nunique()
                else:
                    i_size = mapping["region_b"].nunique()
            coeff = i_size / min_s
    else:
        coeff = np.nan
    return coeff


names_a = []
names_b = []
tf_coefs = []
cre_coefs = []
target_coefs = []
edge_coefs = []
for i in tqdm(range(len(names))):
    name_a = names[i]
    tf_a = tfs[i]
    cr_a = cres[i]
    gn_a = genes[i]
    ed_a = edges[i]
    for j in range(i, len(names)):
        name_b = names[j]
        tf_b = tfs[j]
        cr_b = cres[j]
        gn_b = genes[j]
        ed_b = edges[j]
        names_a.append(name_a)
        names_b.append(name_b)
        tf_coefs.append(set_ocoef(tf_a, tf_b))
        if name_a.startswith('pando') or name_b.startswith('pando'):
            use_overlap = True
        else:
            use_overlap = False
        cre_coefs.append(set_ocoef(cr_a, cr_b, use_overlap=use_overlap))
        target_coefs.append(set_ocoef(gn_a, gn_b))
        edge_coefs.append(set_ocoef(ed_a, ed_b))


# Store as df
sims = pd.DataFrame()
sims['name_a'] = names_a
sims['name_b'] = names_b
sims['tf_oc'] = tf_coefs
sims['cre_oc'] = cre_coefs
sims['target_oc'] = target_coefs
sims['edge_oc'] = edge_coefs

# Write
stats.to_csv(stat_path, index=False)
sims.to_csv(sim_path, index=False)
