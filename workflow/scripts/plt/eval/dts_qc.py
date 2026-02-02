#!/usr/bin/env python
"""
Dataset Summary QC Figure

Generates a summary figure showing QC statistics for all single-cell multiome
datasets. Each row represents a dataset with 4 columns:
1. UMAP with celltype labels
2. Barplot: number of cells per celltype
3. Boxplot: total counts per celltype (hue=omic)
4. Boxplot: number of features per celltype (hue=omic)
"""

import argparse
import mudata as mu
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import yaml


def parse_args():
    parser = argparse.ArgumentParser(description='Generate dataset QC summary figure')
    parser.add_argument('-c', '--config', required=True, help='Path to config.yaml')
    parser.add_argument('-o', '--output', required=True, help='Output PNG path')
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


def main():
    args = parse_args()

    # Load config
    config = load_config(args.config)
    datasets = get_datasets(config)
    n_datasets = len(datasets)

    # Create figure with 4 columns per dataset row
    fig, axes = plt.subplots(
        n_datasets, 4,
        figsize=(10, 2.25 * n_datasets),
        gridspec_kw={'width_ratios': [1, 1, 1, 1]}
    )

    # Ensure axes is 2D even for single row
    if n_datasets == 1:
        axes = axes.reshape(1, -1)

    # Color palette for omics
    omic_palette = {'RNA': 'purple', 'ATAC': '#ffd700'}

    for i, dts_info in enumerate(datasets):
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
            for j in range(4):
                axes[i, j].text(0.5, 0.5, f'Data not found\n{dat}',
                               ha='center', va='center', transform=axes[i, j].transAxes)
                axes[i, j].set_xticks([])
                axes[i, j].set_yticks([])
            continue

        # Sort celltypes for consistent ordering
        celltypes = sorted(df_nc['celltype'].unique())
        df_nc = df_nc.set_index('celltype').loc[celltypes].reset_index()

        # Column 0: UMAP with celltype labels on data
        ax = axes[i, 0]
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

        # Column 1: Barplot of cells per celltype
        ax = axes[i, 1]
        sns.barplot(
            data=df_nc,
            x='size',
            y='celltype',
            order=celltypes,
            ax=ax,
            color='steelblue'
        )
        ax.set_xlabel('# Cells' if i == n_datasets - 1 else '')
        ax.set_ylabel('')
        ax.tick_params(axis='y', labelsize=8)
        if i == 0:
            ax.set_title('Cells per type', fontsize=10)

        # Aggregate QC data: mean per celltype per omic
        df_qc_agg = df_qc.groupby(['celltype', 'omic'], as_index=False).agg({
            'log1p_total_counts': 'mean',
            'log1p_n_genes_by_counts': 'mean'
        })
        df_qc_agg['omic'] = df_qc_agg['omic'].str.upper()

        # Column 2: Boxplot of mean total counts per omic (with stripplot)
        ax = axes[i, 2]
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
        ax.set_xlabel('Omic' if i == n_datasets - 1 else '')
        ax.set_ylabel('log1p(total counts)')
        if i == 0:
            ax.set_title('Total counts', fontsize=10)

        # Column 3: Boxplot of mean features per omic (with stripplot)
        ax = axes[i, 3]
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
        ax.set_xlabel('Omic' if i == n_datasets - 1 else '')
        ax.set_ylabel('log1p(# features)')
        if i == 0:
            ax.set_title('Features detected', fontsize=10)
        del mdata

    # Add shared legend for omics outside the figure
    #from matplotlib.patches import Patch
    #legend_handles = [Patch(facecolor=omic_palette['rna'], label='RNA'),
    #                  Patch(facecolor=omic_palette['atac'], label='ATAC')]
    #fig.legend(handles=legend_handles, loc='center right', title='Omic',
    #           bbox_to_anchor=(1.02, 0.5), fontsize=9, title_fontsize=10)

    plt.tight_layout()
    #plt.subplots_adjust(right=0.92)
    plt.savefig(args.output, dpi=150, bbox_inches='tight')
    plt.close()


if __name__ == '__main__':
    main()
