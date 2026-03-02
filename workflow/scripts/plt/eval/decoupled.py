import pandas as pd
import glob
import sys
import os
import yaml
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats import ttest_1samp
from statsmodels.stats.multitest import multipletests
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils import read_config, savefigs


config = read_config()
class_dict = config['class_names']
db_dict = config['dbs_names']
task_dict = config['task_names']
dts_dict = config['dts_names']
mth_dict = config['method_names']
palette_raw = config['colors']['nets']

# Read
lst_paths = []
for d_dts in ['pitupair', 'lung', 'embryo']:
    lst_paths.extend(glob.glob(f'anl/metrics/*/*/*/hg38.{d_dts}.all.scores.csv'))
df = []
steps = ['pre', 'c2g', 'tfb', 'mdl']
for path in lst_paths:
    spath = path.split('/')
    class_metric = spath[2]
    task_metric = spath[3]
    name_db = spath[4]
    name_dts = spath[5].split('.')[1]
    tmp = pd.read_csv(path)
    msk = tmp['name'].str.startswith('o_')
    tmp = tmp.loc[~msk, :].copy()
    tmp['class'] = class_dict[class_metric]
    tmp['task'] = task_dict[task_metric]
    tmp['db'] = db_dict[name_db]
    tmp['dts'] = name_dts
    df.append(tmp)
cols = ['class', 'task', 'db', 'dts', 'name', 'prc', 'rcl', 'f01']
df = pd.concat(df, axis=0).loc[:, cols].reset_index(drop=True)
df[steps] = df['name'].str.split('.', expand=True, n=4)

# Get mean per dataset; also save per-dataset data for statistical tests
res = []
step_per_dts = {}
for stp in steps:
    mean_step = df.groupby([stp, 'dts', 'class', 'task'], as_index=False)['f01'].mean()
    mean_step = mean_step.groupby([stp, 'dts', 'class'], as_index=False)['f01'].mean()
    mean_step = mean_step.groupby([stp, 'dts',], as_index=False)['f01'].mean()
    step_per_dts[stp] = mean_step.rename(columns={stp: 'name'}).copy()
    mean_step = mean_step.rename(columns={stp: 'name'})
    mean_step['step'] = stp
    res.append(mean_step)
res = pd.concat(res)

# Compute final mean for plotting, sorted by descending mean F0.1
res_mat = (
    res
    .groupby(['name', 'step'], as_index=False)['f01'].mean()
    .pivot(index='name', columns='step', values='f01').loc[:, steps]
)
res_mat = res_mat.loc[res_mat.mean(axis=1).sort_values(ascending=False).index]

# Statistical test: one-sample t-test per (method, step) vs step mean; FDR per step
pval_mat = pd.DataFrame(np.nan, index=res_mat.index, columns=res_mat.columns)
for stp in steps:
    data_stp = step_per_dts[stp]
    step_mean = data_stp['f01'].mean()
    pvals, methods_in_step = [], []
    for mth in res_mat.index:
        g1 = data_stp.loc[data_stp['name'] == mth, 'f01'].values
        if len(g1) < 2:
            pvals.append(np.nan)
        else:
            _, p = ttest_1samp(g1, popmean=step_mean)
            pvals.append(p)
        methods_in_step.append(mth)
    valid = [(i, p) for i, p in enumerate(pvals) if not np.isnan(p)]
    if valid:
        idxs, raw_ps = zip(*valid)
        _, padj, _, _ = multipletests(raw_ps, method='fdr_bh')
        for idx, p_adj in zip(idxs, padj):
            pval_mat.loc[methods_in_step[idx], stp] = p_adj
print(pval_mat)
sig_mat = pval_mat < 0.05

# Build colors for barplot (raw method names as index)
bar_colors = [palette_raw.get(m, '#888888') for m in res_mat.index]
# Rename index to display names for plotting
res_mat.index = [mth_dict.get(m, m) for m in res_mat.index]
sig_mat.index = res_mat.index

# Layout — 3 rows: std barplot (top) | heatmap + right barplot | colorbar (bottom)
n_methods = len(res_mat)
fig = plt.figure(figsize=(6, max(4, n_methods * 0.38 + 2.5)), dpi=150)
gs = gridspec.GridSpec(3, 2, figure=fig,
                       width_ratios=[4, 1.5],
                       height_ratios=[2, n_methods, 2],
                       hspace=0, wspace=0)
ax_top = fig.add_subplot(gs[0, 0])
ax_hm  = fig.add_subplot(gs[1, 0])
ax_bar = fig.add_subplot(gs[1, 1], sharey=ax_hm)
ax_cbar_host = fig.add_subplot(gs[2, 0])
for _ax in [fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[2, 1])]:
    _ax.set_visible(False)
ax_cbar_host.set_axis_off()

# Heatmap
sns.heatmap(
    res_mat,
    cmap='viridis',
    annot=True,
    fmt='.3f',
    ax=ax_hm,
    cbar=False,
    linewidths=0.5,
)
# Asterisk overlay for significant cells
for i, mth in enumerate(res_mat.index):
    for j, stp in enumerate(res_mat.columns):
        if sig_mat.loc[mth, stp] == True:
            ax_hm.text(j + 0.5, i + 0.80, '*',
                       ha='center', va='center',
                       color='white', fontweight='bold')
ax_hm.set_xlabel('')
ax_hm.set_ylabel('')
ax_hm.tick_params(axis='x', labelrotation=0, labelbottom=True, labeltop=False)

# Colorbar — placed in lower portion of ax_cbar_host, ticks at bottom
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import matplotlib.ticker as mticker
vmin = np.nanmin(res_mat.values)
vmax = np.nanmax(res_mat.values)
sm = ScalarMappable(cmap='viridis', norm=Normalize(vmin=vmin, vmax=vmax))
sm.set_array([])
axins = ax_cbar_host.inset_axes([0.275, 0.05, 0.45, 0.4])
cbar = fig.colorbar(sm, cax=axins, orientation='horizontal')
cbar.set_ticks([vmin, (vmin + vmax) / 2, vmax])
cbar.ax.xaxis.set_major_formatter(mticker.FormatStrFormatter('%.2f'))
cbar.ax.tick_params(top=False, labeltop=False, bottom=True, labelbottom=True)
cbar.ax.xaxis.set_label_position('bottom')
cbar.set_label(r'Mean F$_{0.1}$')

# Right barplot: mean F0.1 per method — aligned to heatmap rows via sharey
# seaborn heatmap centers row i at y = i + 0.5
method_means = res_mat.mean(axis=1)
y_positions = np.arange(n_methods) + 0.5
ax_bar.barh(y_positions, method_means.values, height=0.7,
            color=bar_colors, edgecolor='none')
ax_bar.tick_params(axis='y', left=False, right=False, labelleft=False, labelright=False)
ax_bar.set_xlabel(r'Mean F$_{0.1}$')
ax_bar.set_xlim(0, method_means.max() * 1.2)
ax_bar.axvline(method_means.mean(), color='grey', lw=0.8, linestyle='--')
ax_bar.spines[['top', 'right', 'left', 'bottom']].set_visible(False)

# Top std barplot — x aligned to heatmap columns, bars grow top-to-bottom (inverted y)
step_stds = res_mat.std(axis=0)
x_positions = np.arange(len(steps)) + 0.5
ax_top.bar(x_positions, step_stds.values, color='#aaaaaa', width=0.6)
for j, s in enumerate(step_stds.values):
    ax_top.text(x_positions[j], s + step_stds.max() * 0.04, f'{s:.3f}',
                ha='center', va='bottom')
ax_top.set_xlim(0, len(steps))
ax_top.set_ylim(0, step_stds.max() * 1.4)
ax_top.set_ylabel(r'Std F$_{0.1}$')
ax_top.tick_params(axis='x', labelbottom=False, bottom=False)
ax_top.spines[['top', 'right', 'bottom']].set_visible(False)

os.makedirs('plt/eval', exist_ok=True)
savefigs([fig], 'plt/eval/decoupled.pdf')
