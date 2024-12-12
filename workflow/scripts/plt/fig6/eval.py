import marsilea as ma
import marsilea.plotter as mp
from matplotlib.ticker import MaxNLocator
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd
import numpy as np
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils import read_config, savefigs


def summary_steps(df, palette, thr_padj=0.05):
    mat_f01 = df.pivot(index=['metric', 'task', 'db'], columns=['stp', 'name'], values='f01')
    mat_padj = df.pivot(index=['metric', 'task', 'db'], columns=['stp', 'name'], values='padj')
    # Main heatmap
    h1 = ma.Heatmap(
        mat_f01,
        cmap='Purples',
        vmin=0,
        vmax=1,
        label=r'Mean F$\mathrm{_{0.1}}$'
    )
    
    # Add sign asterisks
    mark_sign = mp.MarkerMesh(mat_padj < thr_padj, color="#DB4D6D", label="Sign.")
    h1.add_layer(mark_sign)
    
    # Add labels
    h1.add_left(mp.Labels([i[-1] for i in mat_f01.index], align="right"))
    mth_names = [i[-1] for i in mat_f01.columns]
    h1.add_bottom(mp.Labels([i[-1] for i in mat_f01.columns], align="right"))
    
    # Add barplots
    top_mths = ((mat_padj < thr_padj) * 1).sum(0).reset_index().set_index('name')[0]
    h1.add_top(ma.plotter.Numbers(top_mths, color=[palette[m] for m in top_mths.index], label='Number of\nsign. metrics'))
    top_mtrs = ((mat_padj < thr_padj) * 1).sum(1).reset_index().set_index('db')[0]
    h1.add_right(ma.plotter.Numbers(top_mtrs, color='#3d007d', label='Number of\nsign. methods'))
    
    # Format
    h1.add_legends("right")
    h1.render()
    plt.close()
    fig = h1.figure
    fig.set_dpi(150)

    return fig


def average_ranks(df_list):
    df = pd.concat(df_list)
    df = df.groupby('name')[['rank', 'f01']].mean().sort_values('f01', ascending=False).reset_index()
    df[['pre', 'p2g', 'tfb', 'mdl']] = df['name'].str.split('.', n=4, expand=True)
    df = df.drop(columns='rank').reset_index(names='rank')
    df['fixed'] = [np.unique(n.split('.')).size == 1 for n in df['name']]
    return df


def read_eval(type_metric, task, resource, case, load_originals=False):
    path = 'anl/metrics/{0}/{1}/{2}/{3}.scores.csv'.format(type_metric, task, resource, case)
    df = pd.read_csv(path).sort_values('f01', ascending=False)
    df[['pre', 'p2g', 'tfb', 'mdl']] = df['name'].str.split('.', n=4, expand=True)
    if not load_originals:
        df = df[~df['pre'].str.startswith('o_')]
    df = df.reset_index(drop=True).reset_index(names='rank')
    df['fixed'] = [np.unique(n.split('.')).size == 1 for n in df['name']]
    return df


def ranking(df, title, palette):
    def cat_mat(df, ax, palette):
        cats = list(df['pre'].unique())
        mat = df[['pre', 'p2g', 'tfb', 'mdl']].copy().T.rename({'p2g': 'c2g'})
        for col in mat.columns:
            mat[col] = [cats.index(i) for i in mat[col]]
        mat_palette = {cats.index(k): v for k, v in palette.items()}
        cmap = sns.color_palette([mat_palette[key] for key in sorted(mat_palette)])
        sns.heatmap(mat, cmap=cmap, yticklabels=True, xticklabels='auto', cbar=None, ax=ax)
        ax.set_xlim(-10, df['rank'].max() + 10)
        ax.set_xlabel('Rank')
    
    def plot_ranks(df, ax, palette=None, title='', plot_lg=False):
        ax.plot(df['rank'], df['f01'], c='gray')
        tmp = df[df['fixed']]
        handles = []
        names = []
        for name, rnk, scr in zip(tmp['pre'], tmp['rank'], tmp['f01']):
            if name == 'random':
                ax.axvline(rnk, ls='--', c='black')
                continue
            if name.startswith('o_'):
                linestyle=':'
            else:
                linestyle='-'
            name = name.replace('o_', '')
            linefmt = palette[name]
            markerline, stemlines, baseline = ax.stem(rnk, scr, linefmt=linefmt, orientation='vertical', bottom=-0.15)
            plt.setp(stemlines, linestyle=linestyle, alpha=0.5)
            plt.setp(markerline, markersize=4)
            if name not in names:
                lgn = Line2D([0], [0], color=linefmt, linestyle=linestyle, marker='o', label=name)
                handles.append(lgn)
                names.append(name)
        if plot_lg:
            ax.legend(handles=handles, loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
        ax.invert_yaxis()
        ax.set_ylabel(r'F$\mathrm{_{0.1}}$')
        ax.set_ylim(1.15, -0.15)
        ax.set_yticks([0, 0.5, 1])
        ax.set_title(title)
    fig, axes = plt.subplots(2, 1, figsize=(3, 1.5), dpi=150, sharex=True)
    plot_ranks(df, axes[0], palette=palette)
    axes[0].invert_yaxis()
    cat_mat(df, axes[1], palette=palette)
    axes[1].xaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
    axes[1].set_xticklabels(axes[1].get_xticks().astype(int), rotation=0)
    fig.subplots_adjust(hspace=0.01)
    fig.suptitle(title, fontsize=9)
    return fig


# Read config
config = read_config()
palette = config['colors']['nets']
mthds = list(config['methods'].keys())
baselines = config['baselines']

# Read
df = pd.read_csv(sys.argv[1])

# Sort
df['metric'] = pd.Categorical(df['metric'], categories=['mech', 'pred', 'prior'], ordered=True)
df['task'] = pd.Categorical(df['task'], categories=['tfa', 'prt', 'omics', 'gsets', 'tfm', 'tfp', 'tfb', 'cre', 'c2g'], ordered=True)
df['db'] = pd.Categorical(df['db'], categories=['knocktf', 'gcre', 'cretf', 'gtf', 'hall', 'kegg', 'reac', 'prog', 'hpa', 'tfmdb', 'intact', 'chipatlas', 'remap2022', 'unibind', 'encode', 'gwascatalogue', 'phastcons', 'promoters', 'zhang21', 'eqtlcatalogue'], ordered=True)
df['stp'] = pd.Categorical(df['stp'], categories=['pre', 'p2g', 'tfb', 'mdl'], ordered=True)
df = df.sort_values(['metric', 'task', 'db', 'stp', 'name'])

# Plot
figs = []
figs.append(summary_steps(df, palette, thr_padj=2.2e-22))
metrics = {
    'mech': {
        'prt': ['knocktf'],
        'tfa': ['knocktf'],
    },
    'pred': {
        'omics': ['gtf', 'cretf', 'gcre'],
        'gsets': ['hall', 'kegg', 'reac', 'prog'],
    },
    'prior': {
        'tfm': ['hpa', 'tfmdb'],
        'tfp': ['intact'],
        'tfb': ['chipatlas', 'remap2022', 'unibind'],
        'cre': ['encode', 'phastcons', 'gwascatalogue', 'promoters', 'zhang21'],
        'c2g': ['eqtlcatalogue']
    },
}

case = 'pbmc10k.all'
for mtrc in metrics:
    mtrc_rnks = []
    for task in metrics[mtrc]:
        task_rnks = []
        for db in metrics[mtrc][task]:
            rnk = read_eval(type_metric=mtrc, task=task, resource=db, case=case)
            title = f'{mtrc}|{task}|{db}'
            figs.append(ranking(rnk.dropna(), title, palette))
            task_rnks.append(rnk)
        task_rnks = average_ranks(task_rnks)
        title = f'{mtrc}|{task}'
        mtrc_rnks.append(task_rnks)
        figs.append(ranking(task_rnks.dropna(), title, palette))
    mtrc_rnks = average_ranks(mtrc_rnks)
    title = f'{mtrc}'
    figs.append(ranking(mtrc_rnks.dropna(), title, palette))


# Write
savefigs(figs, sys.argv[2])
