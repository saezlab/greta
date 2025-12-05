import pandas as pd
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils import f_beta_score

path_out = sys.argv[1]

gt = pd.read_csv('dbs/sim/GT_GRN.csv', index_col=0).rename(columns={'regulator': 'source'})

res = []
names = ['celloracle', 'figr', 'pando', 'pearson', 'spearman', 'grnboost', 'random']
for mth in names:
    grn = pd.read_csv(f'dts/sim/seed_1/{mth}.csv')
    
    set_gt = set(gt['source'] + '|' + gt['target'])
    set_grn = set(grn['source'] + '|' + grn['target'])
    
    tp = len(set_gt & set_grn)
    fp = len(set_grn - set_gt)
    fn = len(set_gt - set_grn)

    prc = tp / (tp + fp)
    rcl = tp / (tp + fn)
    f01 = f_beta_score(prc, rcl)

    res.append([mth, prc, rcl, f01])

res = pd.DataFrame(res, columns=['name', 'prc', 'rcl', 'f01'])
res.to_csv(path_out, index=False)
