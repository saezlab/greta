import pandas as pd
import numpy as np
import os
import sys
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils import ocoeff

def read_config(path_config='config/config.yaml'):
    import yaml
    with open(path_config, 'r') as file:
        config = yaml.safe_load(file)
    return config

path_out = sys.argv[1]
config = read_config()

names = ['celloracle', 'figr', 'pando', 'pearson', 'spearman', 'grnboost', 'random']

grns = [pd.read_csv(f'dts/sim/seed_1/{i}.csv') for i in names]
n = len(names)
df = []
for i in range(n):
    for j in range(n):
        s = ocoeff(grns[i], grns[j], on=['source'])
        e = ocoeff(grns[i], grns[j], on=['source', 'target'])
        t = ocoeff(grns[i], grns[j], on=['target'])
        df.append([names[i], names[j], s, e, t])
df = pd.DataFrame(df, columns=['name_a', 'name_b', 'tf_oc', 'edge_oc', 'target_oc'])
mth_names = config['method_names']
df['name_a'] = [mth_names[m] for m in df['name_a']]
df['name_b'] = [mth_names[m] for m in df['name_b']]

df.to_csv(path_out, index=False)
