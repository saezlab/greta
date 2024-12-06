import seaborn as sns
import marsilea as ma
import marsilea.plotter as mp
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils import read_config, savefigs


def fixed_pip(mthds, sts, mat, title):
    res = []
    steps = ['pre', 'p2g', 'tfb', 'mdl']
    for mth in mthds:
        # Extract steps
        msk_mth = (sts['pre'] == mth) & (sts['p2g'] == mth) & (sts['tfb'] == mth) & (sts['mdl'] == mth)
        msk_pre = (sts['pre'] != mth) & (sts['p2g'] == mth) & (sts['tfb'] == mth) & (sts['mdl'] == mth)
        msk_p2g = (sts['pre'] == mth) & (sts['p2g'] != mth) & (sts['tfb'] == mth) & (sts['mdl'] == mth)
        msk_tfb = (sts['pre'] == mth) & (sts['p2g'] == mth) & (sts['tfb'] != mth) & (sts['mdl'] == mth)
        msk_mdl = (sts['pre'] == mth) & (sts['p2g'] == mth) & (sts['tfb'] == mth) & (sts['mdl'] != mth)
        
        # Build df
        df = pd.concat([
            mat.loc[sts[msk_pre].index, sts[msk_mth].index].assign(step=0),
            mat.loc[sts[msk_p2g].index, sts[msk_mth].index].assign(step=1),
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
    
    fig, axes = plt.subplots(len(mthds), 1, sharex=True, sharey=True, figsize=(2, 6.25), dpi=150, tight_layout=True)
    for i, mth in enumerate(mthds):
        df = res[res['mth'] == mth]
        ax = axes[i]
        if i == 0:
            ax.set_title(title)
        sns.boxplot(data=df, x='step', order=steps, y='ocoeff', ax=ax, color=palette[mth], fill=None, fliersize=0)
        sns.stripplot(data=df, x='step', order=steps, y='ocoeff', ax=ax, color=palette[mth])
        ax.set_ylabel(mth)
        ax.set_ylim(-0.05, 1)
        ax.axhline(y=0.5, ls='--', color='black')
    return fig


def sim_mat(mat, sts, palette):
    h1 = ma.Heatmap(mat, cmap='Purples', width=3, height=3, name="h1", vmax=1, label="Overlap\nCoefficient")
    legend=False
    for cat_col in ['mdl', 'tfb', 'p2g', 'pre']:
        cats = sts[cat_col].to_list()
        cat_colors = mp.Colors(cats, palette={k: v for k, v in palette.items() if k in cats}, label=cat_col)
        if cat_col == 'pre':
            legend=True
        h1.add_top(cat_colors, pad=0, size=0.25, legend=legend)
    h1.add_dendrogram("left", show=False)
    h1.add_dendrogram("top")
    h1.add_legends(
        side="right",
        stack_by='col',
        stack_size=2,
    )
    h1.render()
    plt.close()
    hax = h1.get_ax(name='h1')
    border = Rectangle((0, 0), 1, 1, fill=False, ec=".1", lw=2, transform=hax.transAxes)
    hax.add_artist(border)
    fig = h1.figure
    fig.set_dpi(150)
    return fig


# Read config
config = read_config()
palette = config['colors']['nets']
mthds = list(config['methods'].keys())
baselines = config['baselines']

df = pd.read_csv(sys.argv[1])
sts = pd.read_csv(sys.argv[2])

# Remove original runs and baselines
df = df[~(df['name_a'].str.startswith('o_') | df['name_b'].str.startswith('o_'))]
baselines = ['collectri', 'dorothea', 'random', 'scenic']
df = df[~(df['name_a'].str.split('.', expand=True)[0].isin(baselines) | df['name_b'].str.split('.', expand=True)[0].isin(baselines))]

figs = []
for oc in ['tf_oc', 'edge_oc', 'target_oc']:
    mat = df.dropna().pivot(index='name_a', columns='name_b', values=oc).fillna(0)
    mat = mat + mat.T
    np.fill_diagonal(mat.values, 1)
    t_sts = sts.set_index('name').loc[mat.index]
    t_sts[['pre', 'p2g', 'tfb', 'mdl']] = t_sts.reset_index()['name_a'].str.split('.', n=4, expand=True).values
    figs.append(fixed_pip(mthds, t_sts, mat, title=oc))
    figs.append(sim_mat(mat, t_sts, palette))

# Write
savefigs(figs, sys.argv[3])
