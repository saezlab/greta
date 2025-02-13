import pandas as pd
import numpy as np
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils import read_config


def fixed_pip(mthds, sts, mat, title):
    res = []
    steps = ['pre', 'c2g', 'tfb', 'mdl']
    for mth in mthds:
        # Extract steps
        msk_mth = (sts['pre'] == mth) & (sts['c2g'] == mth) & (sts['tfb'] == mth) & (sts['mdl'] == mth)
        msk_pre = (sts['pre'] != mth) & (sts['c2g'] == mth) & (sts['tfb'] == mth) & (sts['mdl'] == mth)
        msk_c2g = (sts['pre'] == mth) & (sts['c2g'] != mth) & (sts['tfb'] == mth) & (sts['mdl'] == mth)
        msk_tfb = (sts['pre'] == mth) & (sts['c2g'] == mth) & (sts['tfb'] != mth) & (sts['mdl'] == mth)
        msk_mdl = (sts['pre'] == mth) & (sts['c2g'] == mth) & (sts['tfb'] == mth) & (sts['mdl'] != mth)
        
        # Build df
        df = pd.concat([
            mat.loc[sts[msk_pre].index, sts[msk_mth].index].assign(step=0),
            mat.loc[sts[msk_c2g].index, sts[msk_mth].index].assign(step=1),
            mat.loc[sts[msk_tfb].index, sts[msk_mth].index].assign(step=2),
            mat.loc[sts[msk_mdl].index, sts[msk_mth].index].assign(step=3),
        ]).reset_index().rename(columns={'{m}.{m}.{m}.{m}'.format(m=mth): 'ocoeff', 'name_a': 'rest'})
        
        # Format df
        df['rest'] = [n.split('.')[i] for n,i in zip(df['rest'], df['step'])]
        df['step'] = [steps[i] for i in df['step']]
        df['mth'] = mth
        df = df[['mth', 'step', 'rest', 'ocoeff']]
        res.append(df)
    res = pd.concat(res)
    res = res.rename(columns={'ocoeff': title})
    return res


# Read
sim = pd.read_csv(sys.argv[1])
sts = pd.read_csv(sys.argv[2])
config = read_config()
mthds = list(config['methods'].keys())

# Remove original runs and baselines
sim = sim[~(sim['name_a'].str.startswith('o_') | sim['name_b'].str.startswith('o_'))]
sim = sim[(sim['name_a'].str.split('.', expand=True)[0].isin(mthds) & sim['name_b'].str.split('.', expand=True)[0].isin(mthds))]

# Find ocoeffs for fixed vs one step change
df = None
for oc in ['tf_oc', 'edge_oc', 'target_oc']:
    mat = sim.dropna().pivot(index='name_a', columns='name_b', values=oc).fillna(0)
    mat = mat + mat.T
    np.fill_diagonal(mat.values, 1)
    t_sts = sts.set_index('name').loc[mat.index].rename(columns={'p2g': 'c2g'})
    t_sts[['pre', 'c2g', 'tfb', 'mdl']] = t_sts.reset_index()['name_a'].str.split('.', n=4, expand=True).values
    if df is None:
        df = fixed_pip(mthds, t_sts, mat, title=oc)
    else:
        df = pd.merge(df, fixed_pip(mthds, t_sts, mat, title=oc), on=['mth', 'step', 'rest'])

# Write
df.to_csv(sys.argv[3], index=False)
