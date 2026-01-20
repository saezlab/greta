import pandas as pd
import scipy.stats as ss
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils import read_config, ocoeff


# Extract dat and case
org, dat, case = os.path.basename(sys.argv[1]).split('.')[:3]

# Read config
config = read_config()
palette = config['colors']['nets']
mthds = list(config['methods'].keys())
is_modular = [config['methods'][m]['modular'] for m in mthds]
m_mthds = [m for m, is_m in zip(mthds, is_modular) if is_m]

# Compute ocoeff and pearson
df = []
for mth in m_mthds:
    ref = pd.read_csv(f'dts/{org}/{dat}/cases/{case}/runs/o_{mth}.o_{mth}.o_{mth}.o_{mth}.grn.csv')
    net = pd.read_csv(f'dts/{org}/{dat}/cases/{case}/runs/{mth}.{mth}.{mth}.{mth}.grn.csv')
    inter = pd.merge(ref, net, on=['source', 'target'], how='inner')
    s, p = ss.pearsonr(inter['score_x'], inter['score_y'])
    df.append([mth, ocoeff(ref, net, on=['source', 'target']), s, p])
df = pd.DataFrame(df, columns=['mth', 'ocoeff', 'stat', 'pval'])
df['padj'] = ss.false_discovery_control(df['pval'])

# Write
df.to_csv(sys.argv[2], index=False)
