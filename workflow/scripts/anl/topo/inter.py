import pandas as pd
import numpy as np
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-g','--paths_grns', required=True, nargs='+')
parser.add_argument('-b','--baselines', required=True, nargs='+')
parser.add_argument('-p','--min_prop', required=True, type=float)
parser.add_argument('-o','--path_out', required=True)
args = parser.parse_args()

grns = []
blns = []
for grn_path in args.paths_grns:
    name = grn_path.split('.')[-3]
    if name.startswith('o_') and (name not in args.baselines):
        grn = pd.read_csv(grn_path).drop_duplicates(['source', 'target'])
        grn['name'] = name.replace('o_', '')
        grns.append(grn)
    elif name in args.baselines:
        grn = pd.read_csv(grn_path).drop_duplicates(['source', 'target']).drop(columns='cre')
        grn['name'] = name
        blns.append(grn)
        
min_n = np.floor(args.min_prop * len(grns))
grns = pd.concat(grns)
blns = pd.concat(blns)
shared = grns.groupby(['source', 'target'], as_index=False).size().sort_values('size', ascending=False)
shared = shared[shared['size'] > min_n]


shared_grn = (
    pd.merge(grns, shared, how='inner', on=['source', 'target'])
    .sort_values(['name', 'source', 'target', 'pval'])
    [['name', 'source', 'target', 'score']]
)
nodes = set(shared_grn['source']) | set(shared_grn['target'])
msk = blns['source'].isin(nodes) & blns['target'].isin(nodes)
blns = blns.loc[msk, :]

shared_grn = pd.concat([
    shared_grn.assign(type='mth'),
    blns.assign(type='bsl')
])

# Write
shared_grn.to_csv(args.path_out, index=False)
