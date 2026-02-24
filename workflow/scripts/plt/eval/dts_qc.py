#!/usr/bin/env python
"""
Dataset Summary QC Figure

Generates a multi-page PDF with one page per dataset. Each page has 4 panels:
1. UMAP with celltype labels (rasterized)
2. Barplot: number of cells per celltype
3. Boxplot: total counts per celltype (hue=omic)
4. Boxplot: number of features per celltype (hue=omic)
"""

import argparse
import mudata as mu
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import yaml


def parse_args():
    parser = argparse.ArgumentParser(description='Generate dataset QC summary figure')
    parser.add_argument('-c', '--config', required=True, help='Path to config.yaml')
    parser.add_argument('-o', '--output', required=True, help='Output PDF path')
    return parser.parse_args()


def load_config(path_config):
    with open(path_config, 'r') as f:
        return yaml.safe_load(f)


def get_datasets(config):
    """Get list of datasets to include, excluding fakepitupair and pitunpair."""
    exclude = {'fakepitupair', 'pitunpair'}
    datasets = []
    for dat, info in config['dts'].items():
        if dat not in exclude:
            datasets.append({
                'dat': dat,
                'org': info['organism'],
                'name': config['dts_names'][dat]
            })
    return datasets


def plot_dataset(pdf, dts_info, omic_palette):
    """Plot a single dataset as one page in the PDF."""
    org = dts_info['org']
    dat = dts_info['dat']
    name = dts_info['name']

    # Define file paths
    mdata_path = f'dts/{org}/{dat}/cases/all/mdata.h5mu'
    qc_path = f'anl/dts/{org}.{dat}.all.qc.csv'
    nc_path = f'anl/dts/{org}.{dat}.all.nc.csv'

    # Load data
    try:
        mdata = mu.read(mdata_path)
        df_qc = pd.read_csv(qc_path)
        df_nc = pd.read_csv(nc_path)
    except FileNotFoundError as e:
        print(f"Warning: Could not load data for {dat}: {e}")
        return

    # Sort celltypes for consistent ordering
    celltypes = sorted(df_nc['celltype'].unique())
    n_celltypes = len(celltypes)
    df_nc = df_nc.set_index('celltype').loc[celltypes].reset_index()

    # Scale height based on number of celltypes (for the barplot)
    bar_height = max(2.5, 0.3 * n_celltypes)
    fig_height = max(3, bar_height)

    fig, axes = plt.subplots(
        1, 4,
        figsize=(12, fig_height),
        gridspec_kw={'width_ratios': [1.2, 1, 0.8, 0.8], 'wspace': 0.4}
    )

    # Panel 0: UMAP with celltype labels (rasterized)
    ax = axes[0]
    sc.pl.embedding(
        mdata,
        basis='X_umap',
        color='celltype',
        ax=ax,
        frameon=False,
        show=False,
        legend_loc='on data',
        legend_fontsize=6,
        legend_fontoutline=2,
        title=name,
    )
    # Rasterize the scatter points for smaller file size
    for child in ax.get_children():
        if hasattr(child, 'set_rasterized'):
            child.set_rasterized(True)

    # Panel 1: Barplot of cells per celltype
    ax = axes[1]
    sns.barplot(
        data=df_nc,
        x='size',
        y='celltype',
        order=celltypes,
        ax=ax,
        color='steelblue'
    )
    ax.set_xlabel('# Cells')
    ax.set_ylabel('')
    ax.set_title('Cells per type', fontsize=10)
    ax.tick_params(axis='y', labelsize=8)

    # Aggregate QC data: mean per celltype per omic
    df_qc_agg = df_qc.groupby(['celltype', 'omic'], as_index=False).agg({
        'log1p_total_counts': 'mean',
        'log1p_n_genes_by_counts': 'mean'
    })
    df_qc_agg['omic'] = df_qc_agg['omic'].str.upper()

    # Panel 2: Boxplot of mean total counts per omic
    ax = axes[2]
    sns.boxplot(
        data=df_qc_agg,
        x='omic',
        y='log1p_total_counts',
        hue='omic',
        ax=ax,
        palette=omic_palette,
        fliersize=0,
        linewidth=0.8,
        fill=False,
        legend=False
    )
    sns.stripplot(
        data=df_qc_agg,
        x='omic',
        y='log1p_total_counts',
        hue='omic',
        ax=ax,
        palette=omic_palette,
        size=4,
        legend=False
    )
    ax.set_xlabel('Omic')
    ax.set_ylabel('log1p(total counts)')
    ax.set_title('Total counts', fontsize=10)

    # Panel 3: Boxplot of mean features per omic
    ax = axes[3]
    sns.boxplot(
        data=df_qc_agg,
        x='omic',
        y='log1p_n_genes_by_counts',
        hue='omic',
        ax=ax,
        palette=omic_palette,
        fliersize=0,
        linewidth=0.8,
        fill=False,
        legend=False
    )
    sns.stripplot(
        data=df_qc_agg,
        x='omic',
        y='log1p_n_genes_by_counts',
        hue='omic',
        ax=ax,
        palette=omic_palette,
        size=4,
        alpha=0.7,
        legend=False
    )
    ax.set_xlabel('Omic')
    ax.set_ylabel('log1p(# features)')
    ax.set_title('Features detected', fontsize=10)

    del mdata

    fig.tight_layout()
    pdf.savefig(fig, dpi=300, bbox_inches='tight')
    plt.close(fig)


def main():
    args = parse_args()

    # Load config
    config = load_config(args.config)
    datasets = get_datasets(config)

    # Color palette for omics
    omic_palette = {'RNA': 'purple', 'ATAC': '#ffd700'}

    with PdfPages(args.output) as pdf:
        for dts_info in datasets:
            print(f"Plotting {dts_info['dat']}...")
            plot_dataset(pdf, dts_info, omic_palette)


if __name__ == '__main__':
    main()
