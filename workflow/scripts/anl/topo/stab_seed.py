import pandas as pd
from glob import glob
from tqdm import tqdm
import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils import (
    ocoeff,
    read_config,
)

config = read_config()
mth_dict = config['method_names']

path_out = sys.argv[1]
lst = path_out.split('/')
org, dts, _, _ = os.path.basename(path_out).split('.')

res = []
for mth in tqdm(mth_dict):
    paths = f'dts/{org}/{dts}/cases/16384_16384_*/runs/o_{mth}.o_{mth}.o_{mth}.o_{mth}.grn.csv'
    grns = []
    for path in glob(paths):
        grn = pd.read_csv(path)
        grns.append(grn)
    n_grns = len(grns)
    assert n_grns == 10
    m_name = mth_dict[mth]
    for i in range(n_grns):
        grn_a = grns[i]
        for j in range(i + 1, n_grns):
            grn_b = grns[j]
            oc_source = ocoeff(df_a=grn_a, df_b=grn_b, on=['source'])
            oc_cre = ocoeff(df_a=grn_a, df_b=grn_b, on=['cre'])
            oc_target = ocoeff(df_a=grn_a, df_b=grn_b, on=['target'])
            oc_edge = ocoeff(df_a=grn_a, df_b=grn_b, on=['source', 'target'])
            res.append([m_name, i, j, oc_source, oc_cre, oc_target, oc_edge])
res = pd.DataFrame(res, columns=['name', 'seed_a', 'seed_b', 'oc_source', 'oc_cre', 'oc_target', 'oc_edge'])

# Write
res.to_csv(path_out, index=False)
