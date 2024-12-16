import matplotlib.patches as patches
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import scipy.stats as ss
import decoupler as dc
import pyranges as pr
import seaborn as sns
import mudata as mu
import pandas as pd
import numpy as np
import argparse
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils import read_config, savefigs


def norm_score(x, axis=None):
    if np.all(x == x[0]):
        return x
    return (x - x.min(axis=axis)) / (x.max(axis=axis) - x.min(axis=axis))


def window(gene, w_size):
    strand = gene.df['Strand'].values[0]
    if strand == '+':
        tss = gene.df['Start'].values[0]
    else:
        tss = gene.df['End'].values[0]
    tss_window = pr.from_dict({
        "Chromosome": gene.Chromosome,
        "Start": tss - w_size,
        "End": tss + w_size,
    })
    return tss_window


def get_links(mthds, baselines, dat, case, target):
    links = []
    for mth in mthds:
        tfb = pd.read_csv(f'dts/{dat}/cases/{case}/runs/{mth}.{mth}.{mth}.tfb.csv').rename(columns={'score': 'tfb_score'})
        if np.isinf(tfb['tfb_score']).any():
            max_finite = tfb['tfb_score'][np.isfinite(tfb['tfb_score'])].max()
            tfb['tfb_score'] = tfb['tfb_score'].replace(np.inf, max_finite)
        tfb['tfb_score'] = norm_score(tfb['tfb_score'].values)
         
        p2g = pd.read_csv(f'dts/{dat}/cases/{case}/runs/{mth}.{mth}.p2g.csv').rename(columns={'score': 'p2g_score'})
        p2g['p2g_score'] = norm_score(p2g['p2g_score'].values)
        p2g = p2g[p2g['gene'] == target]
        link = pd.merge(tfb, p2g, on='cre')[['tf', 'cre', 'gene', 'tfb_score', 'p2g_score']]
        link['mth'] = mth
    
        mdl = pd.read_csv(f'dts/{dat}/cases/{case}/runs/{mth}.{mth}.{mth}.{mth}.mdl.csv').rename(columns={'source': 'tf', 'target': 'gene', 'score': 'mdl_score'})
        mdl['mdl_score'] = norm_score(mdl['mdl_score'].abs().values)
        mdl['mdl_score'] = mdl['mdl_score'].abs().rank(method='average', pct=True)
        link = pd.merge(mdl[['tf', 'gene', 'mdl_score']], link, how='inner')
        links.append(link)
    for mth in baselines:
        mdl = pd.read_csv(f'dts/{dat}/cases/{case}/runs/{mth}.{mth}.{mth}.{mth}.grn.csv').rename(columns={'source': 'tf', 'target': 'gene', 'score': 'mdl_score'})
        mdl['mth'] = mth
        mdl = mdl[mdl['gene'] == target]
        mdl['mdl_score'] = mdl['mdl_score'].abs().rank(method='average', pct=True)
        links.append(mdl)
    links = pd.concat(links).drop(columns=['tfb_score', 'p2g_score'])
    links = links.assign(score=lambda x: x['mdl_score'])
    return links


def get_gannot(path_gannot, target, w_size, atac):
    gannot = pr.read_bed(path_gannot)
    target_gr = gannot[gannot.df['Name'] == target]
    wind = window(target_gr, w_size=w_size)
    x_min = wind.df['Start'].values[0]
    x_max = wind.df['End'].values[0]
    gs_gr = gannot.join(wind).sort()
    
    chromosome = gs_gr.df['Chromosome'].values[0]
    strand = target_gr.df['Strand'].values[0]
    if strand == '+':
        tss = target_gr.df['Start'].values[0]
    else:
        tss = target_gr.df['End'].values[0]
    
    atac_vars = atac.var_names[atac.var_names.str.startswith(chromosome)]
    cres_gr = []
    for cre in atac_vars:
        cre_chr, cre_start, cre_end = cre.split('-')
        if cre_chr == chromosome:
            cres_gr.append([cre_chr, cre_start, cre_end, cre])
    cres_gr = pr.PyRanges(pd.DataFrame(cres_gr, columns=['Chromosome', 'Start', 'End', 'Name']))
    
    return x_min, x_max, gs_gr, tss, chromosome, strand, cres_gr


def mean_data(dat, case):
    mdata = mu.read(f'dts/{dat}/cases/{case}/mdata.h5mu')
    rna = mdata.mod['rna']
    rna.obs = mdata.obs
    atac = mdata.mod['atac']
    atac.obs = mdata.obs
    rna = dc.get_pseudobulk(
        adata=rna,
        sample_col='celltype',
        groups_col=None,
        min_cells=0,
        min_counts=0,
        mode='mean'
    )
    atac = dc.get_pseudobulk(
        adata=atac,
        sample_col='celltype',
        groups_col=None,
        min_cells=0,
        min_counts=0,
        mode='mean'
    )
    return rna, atac


def plot_tfb(tfb_gr, ax):
    for s, e, score in zip(tfb_gr.df['Start'], tfb_gr.df['End'], tfb_gr.df['Score']):
        rect = patches.Rectangle((s, 0), width=e-s, height=score, color='gray')
        ax.add_patch(rect)


def plot_links(links, tf, tss, strand, cmap, ax):
    mthds = links['mth'].sort_values().unique()
    links = links.copy()
    links = links[links['tf'] == tf]
    is_empty = []
    for mth in mthds:
        link = links[links['mth'] == mth]
        if (link.shape[0] == 0):
            is_empty.append(True)
        else:
            is_empty.append(False)
        for i, row in link.iterrows():
            cre, score = row['cre'], row['score']
            if score > 0.00:
                _, cre_start, cre_end = cre.split('-')
                cre_start, cre_end = float(cre_start), float(cre_end)
                if strand == '+':
                    if cre_end < tss:
                        arc_start = cre_end
                        arc_width = tss - cre_end
                    elif cre_start > tss:
                        arc_start = tss
                        arc_width = cre_start - tss
                    elif cre_start < tss < cre_end:
                        arc_start = cre_start
                        arc_width = tss - cre_start
                else:
                    if cre_end < tss:
                        arc_start = cre_end
                        arc_width = tss - cre_end
                    elif cre_start > tss:
                        arc_stat = tss
                        arc_width = cre_start - tss
                    elif cre_start < tss < cre_end:
                        arc_start = tss
                        arc_width = cre_end - tss
                center_x = (arc_start + arc_start + arc_width) // 2
                arc = patches.Arc((center_x, 0), arc_width, score * 2, angle=0, theta1=0, theta2=180, color=cmap[mth], linewidth=1)
                ax.add_patch(arc)
    ax.set_ylim(0, 1.05)
    handles = [plt.Line2D([0], [0], marker='o', color='w', label=mth, markerfacecolor=cmap[mth], markersize=10) 
           for mth in mthds]
    handles = [h for i, h in enumerate(handles) if (not is_empty[i])]
    ax.legend(handles=handles, loc='center left', bbox_to_anchor=(1, 0.5), title='', frameon=False)
    ax.set_ylabel(tf)
    yticks = np.arange(0.25, 1, 0.25)
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticks)
    

def plot_gannot(gs_gr, x_min, x_max, ax):
    def plot_gene(row, i, ax):
        row = row[['Chromosome', 'Start', 'End', 'Strand', 'Name']]
        chrm, strt, end, strnd, name = row
        ax.plot([strt, end], [i, i], lw=1, zorder=0, color='#440a82')
        if strnd == '+':
            s_marker = '>'
        else:
            s_marker = '<'
        ax.scatter([strt, end], [i, i], marker=s_marker, color='#440a82', s=10)
        ax.text(x_max, i, name, ha='left', va='center')
    
    def extract_g_coords(row, size=1000):
        row = row[['Chromosome', 'Start', 'End', 'Strand', 'Name']]
        chrm, strt, end, strnd, name = row
        h_size = size // 2
        if strnd == '+':
            tss_start = strt + h_size
            tss_end = strt - h_size
        else:
            tss_start = end + h_size
            tss_end = end - h_size
        return chrm, tss_start, tss_end, name
        
    
    g_coords = []
    for i, row in gs_gr.df.iterrows():
        plot_gene(row, i, ax)
        g_coords.append(extract_g_coords(row, size=1000))
    g_coords = pr.PyRanges(pd.DataFrame(g_coords, columns=['Chromosome', 'Start', 'End', 'Name']))
    
    ax.set_xlim(x_min, x_max)
    ax.tick_params(axis='y', labelleft=False, length=0)
    w_size = x_max - x_min
    xticks = np.arange(x_min, x_max + 1, w_size // 2)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticks)
    pad =  gs_gr.df.shape[0] * 0.1
    ax.set_ylim(0 - pad, gs_gr.df.shape[0] + pad - 1)


def plot_omic(adata, feat_gr, x_min, x_max, cmap, mode, ax):
    for y, ctype in enumerate(adata.obs_names):
        if mode == 'heatmap':
            y += 0.5
            ax.text(x_max, y, ctype, ha='left', va='center')
            for s, e, feat_name in zip(feat_gr.df['Start'], feat_gr.df['End'], feat_gr.df['Name']):
                if feat_name in adata.var_names:
                    val = adata[ctype, feat_name].X.ravel()[0]
                    rect = patches.Rectangle((s, y - 0.5), width=e-s, height=0.95, color=cmap(val))
                    ax.add_patch(rect)
        elif mode == 'peaks':
            x_coord = [x_min]
            y_coord = [y]
            for s, e, feat_name in zip(feat_gr.df['Start'], feat_gr.df['End'], feat_gr.df['Name']):
                if feat_name in adata.var_names:
                    if s >= x_min and e <= x_max:
                        val = adata[ctype, feat_name].X.ravel()[0]
                        x_coord.extend([s, e])
                        y_coord.extend([val + y, y])
            x_coord.append(x_max)
            y_coord.append(y)
            y_coord = ((norm_score(np.array(y_coord)) * 0.75) + y) + 0.05
            ax.step(x_coord, y_coord, where='post', color='gold')
                
    ax.grid(axis='y', lw=0.25)
    yticks = np.arange(0, adata.shape[0] + 1)
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticks)
    ax.tick_params(axis='y', labelleft=False, length=0)
    ax.set_ylim(0, y + 1)


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-s','--path_sims', required=True)
parser.add_argument('-g','--target', required=True)
parser.add_argument('-t','--tfs', required=True, nargs='+')
parser.add_argument('-a','--path_gannot', required=True)
parser.add_argument('-w','--w_size', required=True, type=int)
parser.add_argument('-o','--path_out', required=True)
args = vars(parser.parse_args())

path_sims = args['path_sims']
target = args['target']
tfs = args['tfs']
path_gannot = args['path_gannot']
w_size = args['w_size']
path_out = args['path_out']

# Extract dataset and case
dat, case = os.path.basename(path_sims).split('.')[:2]

# Read config
config = read_config()
palette = config['colors']['nets']
mthds = list(config['methods'].keys())
baselines = config['baselines']

# Find links
links = get_links(mthds, baselines, dat, case, target)

# Summarize per celltype
rna, atac = mean_data(dat, case)

# Extract gannot data
x_min, x_max, gs_gr, tss, chromosome, strand, cres_gr = get_gannot(path_gannot, target, w_size, atac)

# Plot
fig, axes = plt.subplots(2 + len(tfs), 1, figsize=(3, len(tfs) + 3), dpi=150, sharex=True, height_ratios=[1 for i in range(len(tfs))] + [2, 1])
axes = axes.ravel()
colors = ['white', '#3f007d']
cmap = mcolors.LinearSegmentedColormap.from_list('custom_cmap', colors)
for i in range(len(tfs)):
    ax = axes[i]
    plot_links(links, tfs[i], tss, strand, palette, ax)
ax = axes[-2]
plot_omic(rna, gs_gr, x_min, x_max, cmap, 'heatmap', ax)
plot_omic(atac, cres_gr, x_min, x_max, cmap, 'peaks', ax)
ax = axes[-1]
plot_gannot(gs_gr, x_min, x_max, ax)
ax.set_xlabel(chromosome)
fig.subplots_adjust(wspace=0, hspace=0.0)

# Write
savefigs([fig], path_out)
