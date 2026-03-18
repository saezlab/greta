import pandas as pd
import glob
import sys
import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mticker
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from scipy.stats import ttest_1samp
from statsmodels.stats.multitest import multipletests
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils import read_config, savefigs


config = read_config()
class_dict = config['class_names']
db_dict = config['dbs_names']
task_dict = config['task_names']
mth_dict = config['method_names']
palette_raw = config['colors']['nets']

# Read
lst_paths = []
for d_dts in ['pitupair', 'lung', 'embryo', 'pbmc10k', 'skin']:
    lst_paths.extend(glob.glob(f'anl/metrics/*/*/*/hg38.{d_dts}.all.scores.csv'))
steps = ['pre', 'c2g', 'tfb', 'mdl']
df = []
for path in lst_paths:
    spath = path.split('/')
    tmp = pd.read_csv(path)
    msk = tmp['name'].str.startswith('o_')
    tmp = tmp.loc[~msk, :].copy()
    tmp['class'] = class_dict[spath[2]]
    tmp['task'] = task_dict[spath[3]]
    tmp['db'] = db_dict[spath[4]]
    tmp['dts'] = spath[5].split('.')[1]
    df.append(tmp)
cols = ['class', 'task', 'db', 'dts', 'name', 'prc', 'rcl', 'f01']
df = pd.concat(df, axis=0).loc[:, cols].dropna().reset_index(drop=True)
df[steps] = df['name'].str.split('.', expand=True, n=4)


def compute_matrices(df_sub):
    res = []
    step_per_dts = {}
    for stp in steps:
        mean_step = df_sub.groupby([stp, 'dts', 'class', 'task'], as_index=False)['f01'].mean()
        mean_step = mean_step.groupby([stp, 'dts', 'class'], as_index=False)['f01'].mean()
        mean_step = mean_step.groupby([stp, 'dts'], as_index=False)['f01'].mean()
        step_per_dts[stp] = mean_step.rename(columns={stp: 'name'}).copy()
        mean_step = mean_step.rename(columns={stp: 'name'})
        mean_step['step'] = stp
        res.append(mean_step)
    res_mat = (
        pd.concat(res)
        .groupby(['name', 'step'], as_index=False)['f01'].mean()
        .pivot(index='name', columns='step', values='f01').loc[:, steps]
    )
    res_mat = res_mat.loc[res_mat.mean(axis=1).sort_values(ascending=False).index]
    return res_mat, step_per_dts


def make_fig(df_sub, title):
    res_mat, step_per_dts = compute_matrices(df_sub)

    # One-sample t-test vs leave-one-out mean; FDR per step
    pval_mat = pd.DataFrame(np.nan, index=res_mat.index, columns=res_mat.columns)
    for stp in steps:
        data_stp = step_per_dts[stp]
        pvals, methods_in_step = [], []
        for mth in res_mat.index:
            # 5 observations: per-dataset mean f01 for this method & step
            g1 = data_stp.loc[data_stp['name'] == mth]['f01'].values
            # background: other methods, same step — per-dataset mean then averaged
            other = data_stp.loc[data_stp['name'] != mth]
            step_mean = other.groupby('dts')['f01'].mean().mean()
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
    sig_mat = pval_mat < 0.1

    # Per-figure color scale
    global_vmin = np.nanmin(res_mat.values)
    global_vmax = np.nanmax(res_mat.values)

    # Colors and display names
    bar_colors = [palette_raw.get(m, '#888888') for m in res_mat.index]
    res_mat.index = [mth_dict.get(m, m) for m in res_mat.index]
    sig_mat.index = res_mat.index

    # Layout — 3 rows: std barplot (top) | heatmap + right barplot | colorbar (bottom)
    n_methods = len(res_mat)
    fig = plt.figure(figsize=(7, max(4, n_methods * 0.38 + 2.5)), dpi=150)
    gs = gridspec.GridSpec(3, 2, figure=fig,
                           width_ratios=[4, 2],
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
    sns.heatmap(res_mat, cmap='viridis', annot=True, fmt='.3f',
                ax=ax_hm, cbar=False, linewidths=0.5,
                vmin=global_vmin, vmax=global_vmax)
    for i, mth in enumerate(res_mat.index):
        for j, stp in enumerate(res_mat.columns):
            if sig_mat.loc[mth, stp] == True:
                ax_hm.text(j + 0.5, i + 0.80, '*',
                           ha='center', va='center',
                           color='red', fontweight='bold')
    ax_hm.set_xlabel('')
    ax_hm.set_ylabel('')
    ax_hm.tick_params(axis='x', labelrotation=0, labelbottom=True, labeltop=False)

    # Colorbar
    sm = ScalarMappable(cmap='viridis', norm=Normalize(vmin=global_vmin, vmax=global_vmax))
    sm.set_array([])
    axins = ax_cbar_host.inset_axes([0.275, 0.05, 0.45, 0.4])
    cbar = fig.colorbar(sm, cax=axins, orientation='horizontal')
    cbar.set_ticks([global_vmin, (global_vmin + global_vmax) / 2, global_vmax])
    cbar.ax.xaxis.set_major_formatter(mticker.FormatStrFormatter('%.2f'))
    cbar.ax.tick_params(top=False, labeltop=False, bottom=True, labelbottom=True)
    cbar.ax.xaxis.set_label_position('bottom')
    cbar.set_label(r'Mean F$_{0.1}$')

    # Right barplot
    method_means = res_mat.mean(axis=1)
    y_positions = np.arange(n_methods) + 0.5
    ax_bar.barh(y_positions, method_means.values, height=0.7,
                color=bar_colors, edgecolor='none')
    ax_bar.tick_params(axis='y', left=False, right=False, labelleft=False, labelright=False)
    ax_bar.set_xlabel(r'Mean F$_{0.1}$')
    ax_bar.set_xlim(0, method_means.max() * 1.2)
    ax_bar.axvline(method_means.mean(), color='grey', lw=0.8, linestyle='--')
    ax_bar.spines[['top', 'right', 'left', 'bottom']].set_visible(False)

    # Top std barplot
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
    ax_top.set_title(title)

    return fig


def compute_delta_matrices(df_sub):
    # Identify homogeneous base methods (all 4 steps equal)
    homo_mask = (
        (df_sub['pre'] == df_sub['c2g']) &
        (df_sub['c2g'] == df_sub['tfb']) &
        (df_sub['tfb'] == df_sub['mdl'])
    )
    base_vals = df_sub.loc[homo_mask, 'pre'].unique()

    delta_records = []
    for base_val in base_vals:
        base_data = df_sub[homo_mask & (df_sub['pre'] == base_val)]
        base_g = base_data.groupby(['dts', 'class', 'task'], as_index=False)['f01'].mean()
        for stp in steps:
            other_steps = [s for s in steps if s != stp]
            alt_mask = np.ones(len(df_sub), dtype=bool)
            for ots in other_steps:
                alt_mask &= (df_sub[ots] == base_val).values
            alt_mask &= (df_sub[stp] != base_val).values
            alt_data = df_sub[alt_mask]
            for alt_val in alt_data[stp].unique():
                spec_alt = alt_data[alt_data[stp] == alt_val]
                alt_g = spec_alt.groupby(['dts', 'class', 'task'], as_index=False)['f01'].mean()
                merged = pd.merge(base_g, alt_g, on=['dts', 'class', 'task'], suffixes=('_base', '_alt'))
                if merged.empty:
                    continue
                merged['delta'] = merged['f01_alt'] - merged['f01_base']
                merged['alt'] = alt_val
                merged['step'] = stp
                merged['base'] = base_val
                delta_records.append(merged[['dts', 'class', 'task', 'base', 'alt', 'step', 'delta']])

    if not delta_records:
        return None, None

    delta_df = pd.concat(delta_records, axis=0).reset_index(drop=True)

    # For t-test: per (alt, step, dts), avg over class & task, then over base methods
    step_per_dts = {}
    for stp in steps:
        stp_data = delta_df[delta_df['step'] == stp]
        per_dts = (
            stp_data.groupby(['alt', 'dts', 'base', 'class'], as_index=False)['delta'].mean()
                    .groupby(['alt', 'dts', 'base'], as_index=False)['delta'].mean()
                    .groupby(['alt', 'dts'], as_index=False)['delta'].mean()
        )
        step_per_dts[stp] = per_dts

    # Build heatmap matrix: avg over dts
    res_rows = []
    for stp in steps:
        mean_stp = step_per_dts[stp].groupby('alt', as_index=False)['delta'].mean()
        mean_stp['step'] = stp
        res_rows.append(mean_stp)

    res_mat = (
        pd.concat(res_rows)
        .pivot(index='alt', columns='step', values='delta')
        .reindex(columns=steps)
    )
    res_mat = res_mat.loc[res_mat.mean(axis=1).sort_values(ascending=False).index]
    return res_mat, step_per_dts


def make_delta_fig(df_sub, title):
    res_mat, step_per_dts = compute_delta_matrices(df_sub)
    if res_mat is None:
        return None

    # One-sample t-test vs 0; FDR per step
    pval_mat = pd.DataFrame(np.nan, index=res_mat.index, columns=res_mat.columns)
    for stp in steps:
        data_stp = step_per_dts[stp]
        pvals, alts_in_step = [], []
        for alt in res_mat.index:
            g1 = data_stp.loc[data_stp['alt'] == alt].groupby('dts')['delta'].mean().values
            if len(g1) < 2:
                pvals.append(np.nan)
            else:
                _, p = ttest_1samp(g1, popmean=0)
                pvals.append(p)
            alts_in_step.append(alt)
        valid = [(i, p) for i, p in enumerate(pvals) if not np.isnan(p)]
        if valid:
            idxs, raw_ps = zip(*valid)
            _, padj, _, _ = multipletests(raw_ps, method='fdr_bh')
            for idx, p_adj in zip(idxs, padj):
                pval_mat.loc[alts_in_step[idx], stp] = p_adj
    sig_mat = pval_mat < 0.1

    # Symmetric color scale around 0
    abs_max = max(abs(np.nanmin(res_mat.values)), abs(np.nanmax(res_mat.values)))
    global_vmin, global_vmax = -abs_max, abs_max

    # Colors and display names
    bar_colors = [palette_raw.get(m, '#888888') for m in res_mat.index]
    res_mat.index = [mth_dict.get(m, m) for m in res_mat.index]
    sig_mat.index = res_mat.index

    # Layout — same structure as make_fig
    n_methods = len(res_mat)
    fig = plt.figure(figsize=(7, max(4, n_methods * 0.38 + 2.5)), dpi=150)
    gs = gridspec.GridSpec(3, 2, figure=fig,
                           width_ratios=[4, 2],
                           height_ratios=[2, n_methods, 2],
                           hspace=0, wspace=0)
    ax_top = fig.add_subplot(gs[0, 0])
    ax_hm  = fig.add_subplot(gs[1, 0])
    ax_bar = fig.add_subplot(gs[1, 1], sharey=ax_hm)
    ax_cbar_host = fig.add_subplot(gs[2, 0])
    for _ax in [fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[2, 1])]:
        _ax.set_visible(False)
    ax_cbar_host.set_axis_off()

    # Heatmap — custom annotations to handle potential NaN cells
    annot = res_mat.applymap(lambda x: f'{x:.3f}' if pd.notna(x) else '')
    sns.heatmap(res_mat, cmap='PiYG', annot=annot, fmt='',
                ax=ax_hm, cbar=False, linewidths=0.5,
                vmin=global_vmin, vmax=global_vmax)
    for i, alt in enumerate(res_mat.index):
        for j, stp in enumerate(res_mat.columns):
            if sig_mat.loc[alt, stp] == True:
                ax_hm.text(j + 0.5, i + 0.80, '*',
                           ha='center', va='center',
                           color='red', fontweight='bold')
    ax_hm.set_xlabel('')
    ax_hm.set_ylabel('')
    ax_hm.tick_params(axis='x', labelrotation=0, labelbottom=True, labeltop=False)

    # Colorbar
    sm = ScalarMappable(cmap='PiYG', norm=Normalize(vmin=global_vmin, vmax=global_vmax))
    sm.set_array([])
    axins = ax_cbar_host.inset_axes([0.275, 0.05, 0.45, 0.4])
    cbar = fig.colorbar(sm, cax=axins, orientation='horizontal')
    cbar.set_ticks([global_vmin, 0, global_vmax])
    cbar.ax.xaxis.set_major_formatter(mticker.FormatStrFormatter('%.2f'))
    cbar.ax.tick_params(top=False, labeltop=False, bottom=True, labelbottom=True)
    cbar.ax.xaxis.set_label_position('bottom')
    cbar.set_label(r'Mean $\Delta$F$_{0.1}$')

    # Right barplot — diverging, centered at 0
    method_means = res_mat.mean(axis=1)
    y_positions = np.arange(n_methods) + 0.5
    ax_bar.barh(y_positions, method_means.values, height=0.7,
                color=bar_colors, edgecolor='none')
    ax_bar.tick_params(axis='y', left=False, right=False, labelleft=False, labelright=False)
    ax_bar.set_xlabel(r'Mean $\Delta$F$_{0.1}$')
    max_abs_bar = max(abs(method_means.min()), abs(method_means.max()))
    ax_bar.set_xlim(-max_abs_bar * 1.3, max_abs_bar * 1.3)
    ax_bar.axvline(0, color='grey', lw=0.8, linestyle='--')
    ax_bar.spines[['top', 'right', 'left', 'bottom']].set_visible(False)

    # Top std barplot
    step_stds = res_mat.std(axis=0)
    smax = np.nanmax(step_stds.values) if not np.all(np.isnan(step_stds.values)) else 0.01
    x_positions = np.arange(len(steps)) + 0.5
    ax_top.bar(x_positions, step_stds.values, color='#aaaaaa', width=0.6)
    for j, s in enumerate(step_stds.values):
        if not np.isnan(s):
            ax_top.text(x_positions[j], s + smax * 0.04, f'{s:.3f}',
                        ha='center', va='bottom')
    ax_top.set_xlim(0, len(steps))
    ax_top.set_ylim(0, smax * 1.4)
    ax_top.set_ylabel(r'Std $\Delta$F$_{0.1}$')
    ax_top.tick_params(axis='x', labelbottom=False, bottom=False)
    ax_top.spines[['top', 'right', 'bottom']].set_visible(False)
    ax_top.set_title(title)

    return fig


class_order = ['Genomic', 'Predictive', 'Literature', 'Mechanistic']
subsets = [(cls, df[df['class'] == cls]) for cls in class_order] + [('All', df)]

os.makedirs('plt/eval', exist_ok=True)
figs = [make_fig(df_sub, title) for title, df_sub in subsets]
savefigs(figs, 'plt/eval/decoupled.pdf')
figs_delta = [f for f in (make_delta_fig(df_sub, title) for title, df_sub in subsets) if f is not None]
savefigs(figs_delta, 'plt/eval/decoupled_delta.pdf')
