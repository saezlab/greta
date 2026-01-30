#!/usr/bin/env python
"""
GRN Benchmark Ranking Figure

Generates a comprehensive ranking figure for gene regulatory network inference
benchmark results, including a barplot of overall performance and heatmaps
showing per-class and per-dataset rankings.
"""

import argparse
import textwrap
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec
from scipy.stats import rankdata, kruskal, norm, pearsonr
from itertools import combinations
import yaml


def load_data(path_inp, config_path='config/config.yaml'):
    """Load metrics data and config."""
    # Load config
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)

    # Load metrics
    df = pd.read_csv(path_inp)

    # Exclude control datasets
    df = df[~df['dts'].isin(['Synthetic Pituitary', 'Unpaired Pituitary'])]

    return df, config


def load_scalability_data(path_scalability):
    """Load scalability data from CSV.

    Args:
        path_scalability: Path to scalability CSV with columns: mth, h, gb, use_gpu

    Returns:
        DataFrame indexed by method name
    """
    df = pd.read_csv(path_scalability)
    df = df.set_index('mth')
    return df


def load_pair_data(path_pair):
    """Load paired comparison data from CSV."""
    return pd.read_csv(path_pair)


def load_stats_data(path_stats):
    """Load topology statistics data."""
    return pd.read_csv(path_stats)


def compute_scalability_rankings(scalability_df):
    """Compute rankings for scalability metrics where lower is better.

    Args:
        scalability_df: DataFrame with columns 'h' (time) and 'gb' (memory)

    Returns:
        DataFrame with rankings (rank 1 = fastest/least memory)
    """
    rankings = pd.DataFrame(index=scalability_df.index)
    for col in ['h', 'gb']:
        col_values = scalability_df[col].values
        mask = ~np.isnan(col_values)
        if mask.sum() > 0:
            # Lower values get lower (better) ranks
            valid_ranks = rankdata(col_values[mask], method='average')
            result = np.full(len(col_values), np.nan)
            result[mask] = valid_ranks
            rankings[col] = result
    return rankings


def compute_aggregations(df):
    """Compute mean f01 scores at different aggregation levels."""
    # Overall mean per method
    overall_mean = df.groupby('name')['f01'].mean()

    # Per-class mean
    class_mean = df.groupby(['name', 'class'])['f01'].mean().unstack()

    # Per-dataset mean
    dataset_mean = df.groupby(['name', 'dts'])['f01'].mean().unstack()

    return overall_mean, class_mean, dataset_mean


def compute_rankings(values):
    """
    Compute rankings for each column, where rank 1 = best (highest value).
    Uses average method for ties. NaN values remain as NaN in rankings.
    """
    rankings = pd.DataFrame(index=values.index, columns=values.columns)

    for col in values.columns:
        col_values = values[col].values
        # Create mask for non-NaN values
        mask = ~np.isnan(col_values)
        n_valid = mask.sum()

        if n_valid > 0:
            # Rank only non-NaN values (highest = rank 1)
            valid_ranks = n_valid + 1 - rankdata(col_values[mask], method='average')
            # Place ranks back, keeping NaN where original was NaN
            result = np.full(len(col_values), np.nan)
            result[mask] = valid_ranks
            rankings[col] = result
        else:
            rankings[col] = np.nan

    return rankings


def create_figure(overall_mean, class_mean, dataset_mean, class_ranks, dataset_ranks,
                  scalability_df, scalability_ranks, config, use_rank_barplot=False,
                  mean_class_rank=None):
    """Create the ranking figure with barplot, heatmaps, and scalability panel."""
    # Get method colors and baselines from config
    method_colors_old = config['colors']['nets']
    method_colors = {config['method_names'][k]: method_colors_old[k] for k in method_colors_old}
    baselines = set([config['method_names'][b] for b in config['baselines']])

    # Sort methods by overall mean (descending) or by mean rank (ascending) if flag is set
    if use_rank_barplot and mean_class_rank is not None:
        # Sort by mean rank ascending (lowest rank = best = top)
        method_order = mean_class_rank.sort_values(ascending=True).index.tolist()
    else:
        method_order = overall_mean.sort_values(ascending=False).index.tolist()

    # Sort datasets by mean performance across methods (descending)
    dataset_order = dataset_mean.mean().sort_values(ascending=False).index.tolist()

    # Sort classes by mean performance across methods (descending)
    class_order = class_mean.mean().sort_values(ascending=False).index.tolist()

    # Reorder dataframes
    overall_mean = overall_mean.loc[method_order]
    class_mean = class_mean.loc[method_order, class_order]
    dataset_mean = dataset_mean.loc[method_order, dataset_order]
    class_ranks = class_ranks.loc[method_order, class_order]
    dataset_ranks = dataset_ranks.loc[method_order, dataset_order]

    # Reorder scalability data to match method order
    scalability_df = scalability_df.loc[method_order]
    scalability_ranks = scalability_ranks.loc[method_order]

    # Create figure with constrained layout
    fig = plt.figure(figsize=(10, 6), constrained_layout=True) # 14, 6

    # GridSpec layout with 2 rows:
    # Row 1: [names, barplot, gap, class_heatmap, gap, dataset_heatmap, gap, scal_time, scal_mem, scal_gpu]
    # Row 2: [empty spaces and colorbars for f01, time, memory]
    # Ratios based on user request: barplot=2, class=1, dataset=3, scalability=1 (split into 3)
    gs = GridSpec(2, 11, figure=fig,
                  width_ratios=[1, 1, 0.01, 1, 0.01, 3, 0.01, 0.33, 0.33, 0.33, 0.34],
                  height_ratios=[20, 1],
                  wspace=0.02, hspace=0.15)

    ax_names = fig.add_subplot(gs[0, 0])
    ax_bar = fig.add_subplot(gs[0, 1])
    ax_class = fig.add_subplot(gs[0, 3])
    ax_dataset = fig.add_subplot(gs[0, 5])
    ax_scalability = fig.add_subplot(gs[0, 7:11])

    # Colorbar axes in bottom row - balanced sizes
    # F0.1 under class (ratio 1), stab/time/mem under their columns (ratio 0.33 each)
    ax_cbar_f01 = fig.add_subplot(gs[1, 3])
    ax_cbar_stab = fig.add_subplot(gs[1, 7])
    ax_cbar_time = fig.add_subplot(gs[1, 8])
    ax_cbar_mem = fig.add_subplot(gs[1, 9])

    n_methods = len(method_order)
    y_positions = np.arange(n_methods)

    # --- Method names (ax_names) ---
    ax_names.set_xlim(0, 1)
    ax_names.set_ylim(-0.5, n_methods - 0.5)
    ax_names.invert_yaxis()
    ax_names.axis('off')

    for i, method in enumerate(method_order):
        fontweight = 'normal' if method in baselines else 'bold'
        ax_names.text(1, i, method, ha='right', va='center',
                     fontsize=10, color='black', fontweight=fontweight)

    # --- Barplot (ax_bar) ---
    bar_colors = [method_colors.get(m, 'gray') for m in method_order]

    if use_rank_barplot and mean_class_rank is not None:
        # Invert ranks so longer bars = better (transform: max_rank - rank + 1)
        mean_class_rank_ordered = mean_class_rank.loc[method_order]
        display_values = n_methods + 1 - mean_class_rank_ordered.values
        ax_bar.barh(y_positions, display_values, color=bar_colors, height=0.95)
        ax_bar.set_ylim(-0.5, n_methods - 0.5)
        ax_bar.invert_yaxis()
        ax_bar.set_xlabel('Inverted Mean Rank', fontsize=10)
        ax_bar.set_yticks([])
        ax_bar.set_xlim(0, n_methods + 0.5)
    else:
        ax_bar.barh(y_positions, overall_mean.values, color=bar_colors, height=0.95)
        ax_bar.set_ylim(-0.5, n_methods - 0.5)
        ax_bar.invert_yaxis()
        ax_bar.set_xlabel('Mean F0.1', fontsize=10)
        ax_bar.set_yticks([])
        ax_bar.set_xlim(0, overall_mean.max() * 1.05)

    ax_bar.spines['top'].set_visible(False)
    ax_bar.spines['right'].set_visible(False)
    ax_bar.spines['left'].set_visible(False)
    ax_bar.axvline(x=np.mean(display_values), ls='--', lw=1, c='black')

    # --- Shared colormap settings for F0.1 heatmaps ---
    vmin = 0
    vmax = max(class_mean.values.max(), dataset_mean.values.max())
    cmap = plt.cm.viridis

    # Threshold for text color (midpoint)
    threshold = (vmin + vmax) / 2

    # --- Class Heatmap (ax_class) ---
    class_data = class_mean.values
    im_class = ax_class.imshow(class_data, aspect='auto', cmap=cmap, vmin=vmin, vmax=vmax)

    # Add ranking text
    for i in range(n_methods):
        for j in range(len(class_order)):
            rank_val = class_ranks.iloc[i, j]
            cell_val = class_data[i, j]
            # Handle NaN values
            if pd.isna(rank_val) or pd.isna(cell_val):
                rank_str = '-'
                text_color = 'gray'
            else:
                text_color = 'white' if cell_val > threshold else 'black'
                # Format rank as integer if whole number, else show one decimal
                if rank_val == int(rank_val):
                    rank_str = f'{int(rank_val)}'
                else:
                    rank_str = f'{rank_val:.1f}'
            ax_class.text(j, i, rank_str, ha='center', va='center',
                         fontsize=10, color=text_color)

    ax_class.set_xticks(range(len(class_order)))
    class_labels = [config['class_names'].get(c, c) for c in class_order]
    ax_class.set_xticklabels(class_labels, fontsize=10, rotation=90)
    ax_class.xaxis.set_ticks_position('top')
    ax_class.set_yticks([])
    ax_class.set_ylabel('')

    # --- Dataset Heatmap (ax_dataset) ---
    dataset_data = dataset_mean.values
    im_dataset = ax_dataset.imshow(dataset_data, aspect='auto', cmap=cmap, vmin=vmin, vmax=vmax)

    # Add ranking text
    for i in range(n_methods):
        for j in range(len(dataset_order)):
            rank_val = dataset_ranks.iloc[i, j]
            cell_val = dataset_data[i, j]
            # Handle NaN values
            if pd.isna(rank_val) or pd.isna(cell_val):
                rank_str = '-'
                text_color = 'gray'
            else:
                text_color = 'white' if cell_val > threshold else 'black'
                if rank_val == int(rank_val):
                    rank_str = f'{int(rank_val)}'
                else:
                    rank_str = f'{rank_val:.1f}'
            ax_dataset.text(j, i, rank_str, ha='center', va='center',
                           fontsize=10, color=text_color)

    ax_dataset.set_xticks(range(len(dataset_order)))
    ax_dataset.set_xticklabels(dataset_order, fontsize=10, rotation=90)
    ax_dataset.xaxis.set_ticks_position('top')
    ax_dataset.set_yticks([])

    # --- Scalability Heatmap (ax_scalability) ---
    # Prepare scalability data: Stability, Time (h), Mem (GB), GPU
    stab_values = scalability_df['stability'].values
    time_values = scalability_df['h'].values
    mem_values = scalability_df['gb'].values
    gpu_values = scalability_df['use_gpu'].values

    time_ranks = scalability_ranks['h'].values
    mem_ranks = scalability_ranks['gb'].values

    # Normalize for colormaps
    stab_norm = (stab_values - np.nanmin(stab_values)) / (np.nanmax(stab_values) - np.nanmin(stab_values))
    time_norm = (time_values - np.nanmin(time_values)) / (np.nanmax(time_values) - np.nanmin(time_values))
    mem_norm = (mem_values - np.nanmin(mem_values)) / (np.nanmax(mem_values) - np.nanmin(mem_values))

    # Plot stability column with Greens colormap (darker = higher/better)
    stab_data = stab_norm.reshape(-1, 1)
    ax_scalability.imshow(stab_data, aspect='auto', cmap=plt.cm.Greens,
                          vmin=0, vmax=1, extent=[-0.5, 0.5, n_methods - 0.5, -0.5])

    # Plot time column with Reds colormap
    time_data = time_norm.reshape(-1, 1)
    im_time = ax_scalability.imshow(time_data, aspect='auto', cmap=plt.cm.Reds_r,
                                     vmin=0, vmax=1, extent=[0.5, 1.5, n_methods - 0.5, -0.5])

    # Plot memory column with Blues colormap
    mem_data = mem_norm.reshape(-1, 1)
    ax_scalability.imshow(mem_data, aspect='auto', cmap=plt.cm.Blues_r,
                          vmin=0, vmax=1, extent=[1.5, 2.5, n_methods - 0.5, -0.5])

    # Plot GPU column with gray background
    gpu_data = np.zeros((n_methods, 1))
    ax_scalability.imshow(gpu_data, aspect='auto', cmap='Greys', vmin=0, vmax=1,
                          extent=[2.5, 3.5, n_methods - 0.5, -0.5], alpha=0.2)

    # Add text for stability values (2 decimal places)
    for i in range(n_methods):
        stab_val = stab_values[i]
        cell_norm = stab_norm[i]
        if pd.isna(stab_val):
            val_str = '-'
            text_color = 'gray'
        else:
            text_color = 'white' if cell_norm > 0.5 else 'black'
            val_str = f'{stab_val:.2f}'
        ax_scalability.text(0, i, val_str, ha='center', va='center',
                           fontsize=8, color=text_color)

    # Add text for time values (actual hours)
    for i in range(n_methods):
        time_val = time_values[i]
        cell_norm = time_norm[i]
        if pd.isna(time_val):
            val_str = '-'
            text_color = 'gray'
        else:
            text_color = 'white' if cell_norm < 0.5 else 'black'
            # Format: show 1 decimal for values >= 1, 2 decimals for smaller
            if time_val >= 1:
                val_str = f'{time_val:.1f}'
            else:
                val_str = f'{time_val:.2f}'
        ax_scalability.text(1, i, val_str, ha='center', va='center',
                           fontsize=8, color=text_color)

    # Add text for memory values (actual GB)
    for i in range(n_methods):
        mem_val = mem_values[i]
        cell_norm = mem_norm[i]
        if pd.isna(mem_val):
            val_str = '-'
            text_color = 'gray'
        else:
            text_color = 'white' if cell_norm < 0.5 else 'black'
            # Format: show integer if >= 10, 1 decimal otherwise
            if mem_val >= 10:
                val_str = f'{mem_val:.0f}'
            else:
                val_str = f'{mem_val:.1f}'
        ax_scalability.text(2, i, val_str, ha='center', va='center',
                           fontsize=8, color=text_color)

    # Add GPU symbols (checkmark or X) - bold and larger
    for i in range(n_methods):
        gpu_val = gpu_values[i]
        symbol = '✓' if gpu_val else ''
        text_color = 'green' if gpu_val else 'gray'
        ax_scalability.text(3, i, symbol, ha='center', va='center',
                           fontsize=14, fontweight='bold', color=text_color)

    # Configure scalability axes
    ax_scalability.set_xlim(-0.5, 3.5)
    ax_scalability.set_ylim(n_methods - 0.5, -0.5)
    ax_scalability.set_xticks([0, 1, 2, 3])
    ax_scalability.set_xticklabels(['Stability', 'Time (h)', 'Mem (GB)', 'GPU'], fontsize=10, rotation=90)
    ax_scalability.xaxis.set_ticks_position('top')
    ax_scalability.set_yticks([])

    # --- Colorbars in bottom row ---
    # F0.1 colorbar (for class and dataset heatmaps)
    cbar_f01 = plt.colorbar(im_dataset, cax=ax_cbar_f01, orientation='horizontal')
    cbar_f01.set_label('Mean F0.1', fontsize=9)
    cbar_f01.ax.tick_params(labelsize=8)

    # Stability colorbar (Greens)
    sm_stab = plt.cm.ScalarMappable(cmap=plt.cm.Greens,
                                     norm=plt.Normalize(vmin=np.nanmin(stab_values),
                                                        vmax=np.nanmax(stab_values)))
    cbar_stab = plt.colorbar(sm_stab, cax=ax_cbar_stab, orientation='horizontal')
    cbar_stab.set_label('Mean\nOverlap\nCoefficient', fontsize=8)
    cbar_stab.ax.tick_params(labelsize=7)

    # Time colorbar (Reds)
    sm_time = plt.cm.ScalarMappable(cmap=plt.cm.Reds_r,
                                     norm=plt.Normalize(vmin=np.nanmin(time_values),
                                                        vmax=np.nanmax(time_values)))
    cbar_time = plt.colorbar(sm_time, cax=ax_cbar_time, orientation='horizontal')
    cbar_time.set_label('Time (h)', fontsize=8)
    cbar_time.ax.tick_params(labelsize=7)

    # Memory colorbar (Blues)
    sm_mem = plt.cm.ScalarMappable(cmap=plt.cm.Blues_r,
                                    norm=plt.Normalize(vmin=np.nanmin(mem_values),
                                                       vmax=np.nanmax(mem_values)))
    cbar_mem = plt.colorbar(sm_mem, cax=ax_cbar_mem, orientation='horizontal')
    cbar_mem.set_label('Mem (GB)', fontsize=8)
    cbar_mem.ax.tick_params(labelsize=7)

    return fig


def dunns_test(groups):
    """
    Perform Dunn's test for pairwise comparisons after Kruskal-Wallis.

    Args:
        groups: List of arrays, one per group

    Returns:
        DataFrame with pairwise comparisons and raw p-values
    """
    # Pool all data and compute ranks
    all_data = np.concatenate(groups)
    n_total = len(all_data)
    ranks = rankdata(all_data, method='average')

    # Split ranks back into groups
    group_ranks = []
    idx = 0
    for g in groups:
        group_ranks.append(ranks[idx:idx + len(g)])
        idx += len(g)

    # Compute mean ranks and group sizes
    n_groups = len(groups)
    mean_ranks = [np.mean(r) for r in group_ranks]
    group_sizes = [len(g) for g in groups]

    # Compute tie correction
    unique, counts = np.unique(ranks, return_counts=True)
    tie_sum = np.sum(counts ** 3 - counts)
    tie_correction = 1 - tie_sum / (n_total ** 3 - n_total)

    # Variance of rank sums
    variance = (n_total * (n_total + 1) / 12) * tie_correction

    # Compute z-statistics for all pairs
    results = []
    for i, j in combinations(range(n_groups), 2):
        se = np.sqrt(variance * (1 / group_sizes[i] + 1 / group_sizes[j]))
        z = (mean_ranks[i] - mean_ranks[j]) / se
        p_value = 2 * (1 - norm.cdf(abs(z)))  # Two-tailed
        results.append({
            'group1': i,
            'group2': j,
            'z_stat': z,
            'p_value': p_value
        })

    return pd.DataFrame(results)


def benjamini_hochberg(p_values):
    """
    Apply Benjamini-Hochberg correction for multiple testing.

    Args:
        p_values: Array of p-values

    Returns:
        Array of adjusted p-values
    """
    n = len(p_values)
    sorted_idx = np.argsort(p_values)
    sorted_pvals = np.array(p_values)[sorted_idx]

    # BH adjustment
    adjusted = np.zeros(n)
    for i in range(n):
        rank = i + 1
        adjusted[sorted_idx[i]] = sorted_pvals[i] * n / rank

    # Ensure monotonicity (cumulative minimum from the end)
    adjusted_monotonic = np.zeros(n)
    adjusted_monotonic[sorted_idx[-1]] = min(adjusted[sorted_idx[-1]], 1.0)
    for i in range(n - 2, -1, -1):
        adjusted_monotonic[sorted_idx[i]] = min(
            adjusted[sorted_idx[i]],
            adjusted_monotonic[sorted_idx[i + 1]]
        )
        adjusted_monotonic[sorted_idx[i]] = min(adjusted_monotonic[sorted_idx[i]], 1.0)

    return adjusted_monotonic


def compute_dataset_comparison_stats(df):
    """
    Compute statistical comparison of datasets.

    Args:
        df: DataFrame with columns 'name', 'dts', 'f01'

    Returns:
        kw_stat: Kruskal-Wallis H statistic
        kw_pval: Kruskal-Wallis p-value
        posthoc_df: DataFrame with pairwise comparisons
    """
    # Compute mean F0.1 per method for each dataset
    method_dataset_mean = df.groupby(['name', 'dts'])['f01'].mean().reset_index()

    # Get unique datasets
    datasets = method_dataset_mean['dts'].unique()

    # Create groups for Kruskal-Wallis
    groups = []
    dataset_names = []
    for dts in datasets:
        group_data = method_dataset_mean[method_dataset_mean['dts'] == dts]['f01'].values
        groups.append(group_data)
        dataset_names.append(dts)

    # Kruskal-Wallis test
    kw_stat, kw_pval = kruskal(*groups)

    # Dunn's post-hoc test
    posthoc_raw = dunns_test(groups)

    # Apply Benjamini-Hochberg correction
    posthoc_raw['p_adjusted'] = benjamini_hochberg(posthoc_raw['p_value'].values)

    # Add dataset names
    posthoc_raw['dataset1'] = posthoc_raw['group1'].map(lambda x: dataset_names[x])
    posthoc_raw['dataset2'] = posthoc_raw['group2'].map(lambda x: dataset_names[x])

    return kw_stat, kw_pval, posthoc_raw


def get_significance_symbol(p_value):
    """Return asterisk symbol based on p-value."""
    if p_value < 0.001:
        return '***'
    elif p_value < 0.01:
        return '**'
    elif p_value < 0.05:
        return '*'
    return ''


def create_dataset_boxplot_figure(dataset_mean, posthoc_df, config, dataset_order):
    """
    Create boxplot figure comparing F0.1 across datasets with significance bars.

    Args:
        dataset_mean: DataFrame with mean F0.1 per method for each dataset
        posthoc_df: DataFrame with pairwise comparisons
        config: Configuration dictionary
        dataset_order: List of dataset names in desired order

    Returns:
        matplotlib Figure
    """
    # Prepare data for seaborn boxplot (long format)
    plot_data = []
    for dts in dataset_order:
        values = dataset_mean[dts].dropna().values
        for v in values:
            plot_data.append({'Dataset': dts, 'Mean F0.1': v})
    plot_df = pd.DataFrame(plot_data)

    # Create mapping from dataset name to x position (0-indexed for seaborn)
    dts_to_pos = {dts: i + 1 for i, dts in enumerate(dataset_order)}

    # Create figure
    fig, ax = plt.subplots(figsize=(5, 3))

    # Boxplot using seaborn
    sns.boxplot(data=plot_df, x='Dataset', y='Mean F0.1', ax=ax, width=0.6, order=dataset_order, fill=False)

    # Rotate x-tick labels
    ax.tick_params(axis='x', labelsize=10, rotation=90)
    ax.set_xlabel('')
    ax.set_ylabel('Mean F0.1', fontsize=12)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Get significant comparisons (FDR < 0.05)
    significant = posthoc_df[posthoc_df['p_adjusted'] < 0.05].copy()
    print(significant)

    if len(significant) > 0:
        # Get y-axis limits and compute bar positions
        y_max = max([np.max(d) for d in boxplot_data])
        y_range = y_max - ax.get_ylim()[0]
        bar_height = y_range * 0.03
        bar_gap = y_range * 0.04

        # Sort by distance between groups (draw shorter bars first)
        significant['distance'] = abs(
            significant['dataset1'].map(dts_to_pos) -
            significant['dataset2'].map(dts_to_pos)
        )
        significant = significant.sort_values('distance')

        # Track occupied y levels for each x range
        current_y = y_max + bar_gap

        for _, row in significant.iterrows():
            x1 = dts_to_pos[row['dataset1']]
            x2 = dts_to_pos[row['dataset2']]
            if x1 > x2:
                x1, x2 = x2, x1

            symbol = get_significance_symbol(row['p_adjusted'])

            # Draw the bar
            bar_y = current_y
            ax.plot([x1, x1, x2, x2], [bar_y - bar_height/2, bar_y, bar_y, bar_y - bar_height/2],
                    color='black', linewidth=1)

            # Add asterisks
            ax.text((x1 + x2) / 2, bar_y, symbol, ha='center', va='bottom', fontsize=12)

            current_y += bar_gap + bar_height

        # Adjust y limit to accommodate bars
        ax.set_ylim(ax.get_ylim()[0], current_y + bar_gap)

    plt.tight_layout()
    return fig


def compute_database_aggregation(df):
    """Compute mean f01 scores per (database, task) combination, averaged across all datasets.

    Returns:
        db_mean: DataFrame with methods as index, (db, task) tuples as columns
        db_hierarchy: DataFrame with columns [db, task, class] for grouping
    """
    # Mean F0.1 per method per (database, task) combination (averaged across datasets)
    # This handles cases where the same database appears in multiple tasks
    db_task_mean = df.groupby(['name', 'db', 'task'])['f01'].mean().reset_index()
    db_mean = db_task_mean.pivot(index='name', columns=['db', 'task'], values='f01')

    # Get hierarchy info (class) for each (db, task) combination
    db_hierarchy = df[['db', 'task', 'class']].drop_duplicates()

    return db_mean, db_hierarchy


def get_hierarchical_db_order(db_mean, db_hierarchy, class_order):
    """Order (db, task) combinations by class, then task (by mean perf), then database (by mean perf).

    Args:
        db_mean: DataFrame with (db, task) tuples as columns
        db_hierarchy: DataFrame with db, task, class columns
        class_order: List of class names in desired order ['genom', 'pred', 'prior', 'mech']

    Returns:
        ordered_cols: List of (db, task) tuples in hierarchical order
    """
    # Convert MultiIndex columns to set of tuples for reliable membership checking
    available_cols = set(db_mean.columns.tolist())

    ordered_cols = []
    for cls in class_order:
        # Get (db, task) combinations in this class
        cls_rows = db_hierarchy[db_hierarchy['class'] == cls]

        if cls_rows.empty:
            continue

        # Get tasks in this class, sorted by mean performance
        tasks_in_cls = cls_rows['task'].unique()
        task_means = {}
        for t in tasks_in_cls:
            task_dbs = cls_rows[cls_rows['task'] == t]['db'].tolist()
            task_cols = [(d, t) for d in task_dbs if (d, t) in available_cols]
            if task_cols:
                task_means[t] = db_mean[task_cols].mean().mean()
        sorted_tasks = sorted(tasks_in_cls, key=lambda t: task_means.get(t, 0), reverse=True)

        for task in sorted_tasks:
            # Get databases in this task, sorted by mean performance
            task_dbs = cls_rows[cls_rows['task'] == task]['db'].tolist()
            task_cols = [(d, task) for d in task_dbs if (d, task) in available_cols]
            task_cols_sorted = sorted(task_cols, key=lambda c: db_mean[c].mean(), reverse=True)
            ordered_cols.extend(task_cols_sorted)

    return ordered_cols


def wrap_text_by_words(text, max_width=12):
    """Wrap text at word boundaries only, no mid-word breaks."""
    wrapper = textwrap.TextWrapper(
        width=max_width,
        break_long_words=False,
        break_on_hyphens=False
    )
    lines = wrapper.wrap(text)
    wrapped_text = '\n'.join(lines)
    max_line_width = max(len(line) for line in lines) if lines else 0
    num_lines = len(lines)
    return wrapped_text, max_line_width, num_lines


def precompute_wrapped_titles(task_boundaries_per_class, task_names_map, max_width=12):
    """Pre-compute all wrapped task titles and their dimensions."""
    wrapped_titles = {}
    class_max_line_widths = {}
    global_max_lines = 1

    for cls, boundaries in task_boundaries_per_class.items():
        wrapped_titles[cls] = {}
        max_width_in_class = 0
        max_lines_in_class = 1

        for start, end, task in boundaries:
            label = task_names_map.get(task, task)
            wrapped_text, max_line_width, num_lines = wrap_text_by_words(label, max_width)
            wrapped_titles[cls][task] = wrapped_text
            max_width_in_class = max(max_width_in_class, max_line_width)
            max_lines_in_class = max(max_lines_in_class, num_lines)

        class_max_line_widths[cls] = max_width_in_class
        global_max_lines = max(global_max_lines, max_lines_in_class)

    return wrapped_titles, class_max_line_widths, global_max_lines


def create_database_heatmap_figure(db_mean, db_ranks, db_hierarchy, method_order,
                                    class_order, config):
    """Create heatmap with hierarchical headers (Class -> Task -> Database).

    Layout:
    - Rows: methods (same order as main figure)
    - Columns: databases grouped by task, grouped by class (separate heatmaps per TASK)
    - Top headers: Class labels (spanning), Task labels (per task)
    - Bottom: Database names (rotated 90 degrees)
    """
    baselines = set([config['method_names'][b] for b in config['baselines']])
    n_methods = len(method_order)

    # Columns are (db, task) tuples
    col_list = db_mean.columns.tolist()

    # Get config mappings
    task_names_map = config.get('task_names', {})
    class_names_map = config.get('class_names', {})
    db_names_map = config.get('dbs_names', {})
    class_colors = {'genom': '#e6f3ff', 'pred': '#fff3e6', 'prior': '#e6ffe6', 'mech': '#ffe6f3'}

    # Build task structure: list of (class, task, [(db, task), ...]) in order
    task_list = []  # [(cls, task, col_list), ...]
    for cls in class_order:
        cls_rows = db_hierarchy[db_hierarchy['class'] == cls]
        cls_col_list = [(d, t) for d, t in col_list
                        if any((cls_rows['db'] == d) & (cls_rows['task'] == t))]
        if not cls_col_list:
            continue
        # Group by task within this class
        current_task = None
        task_cols = []
        for d, t in cls_col_list:
            if t != current_task:
                if current_task is not None:
                    task_list.append((cls, current_task, task_cols))
                current_task = t
                task_cols = []
            task_cols.append((d, t))
        if current_task is not None:
            task_list.append((cls, current_task, task_cols))

    # Pre-compute wrapped titles and dimensions for each task
    task_wrapped_titles = {}  # {task: wrapped_text}
    task_title_widths = {}  # {task: max_line_width}
    global_max_lines = 1
    for cls, task, cols in task_list:
        label = task_names_map.get(task, task)
        wrapped_text, max_line_width, num_lines = wrap_text_by_words(label, max_width=12)
        task_wrapped_titles[task] = wrapped_text
        task_title_widths[task] = max_line_width
        global_max_lines = max(global_max_lines, num_lines)

    # Calculate width for each task: max(n_cols, title_width_in_cols)
    # where title_width_in_cols approximates how many columns the title needs
    CHARS_PER_COLUMN = 2.5  # approximate characters that fit per heatmap column
    task_widths = []
    for cls, task, cols in task_list:
        n_cols = len(cols)
        title_chars = task_title_widths[task]
        title_width_in_cols = title_chars / CHARS_PER_COLUMN
        effective_width = max(n_cols, title_width_in_cols)
        task_widths.append(effective_width)

    total_cols = sum(len(cols) for _, _, cols in task_list)

    # Build width_ratios: [method_names, task1, task2, gap, task3, ...]
    # Gaps only between classes
    gap_width = 0.3
    width_ratios = [1.5]  # method names column
    gs_col_mapping = []  # maps task index to gs column
    current_class = None
    for i, (cls, task, cols) in enumerate(task_list):
        if current_class is not None and cls != current_class:
            # Add gap between classes
            width_ratios.append(gap_width)
        width_ratios.append(task_widths[i])
        gs_col_mapping.append(len(width_ratios) - 1)
        current_class = cls

    n_width_cols = len(width_ratios)

    # Dynamic task header height based on max lines
    task_header_height = 0.4 * global_max_lines

    fig = plt.figure(figsize=(max(14, total_cols * 0.5), max(6, n_methods * 0.35)))

    # GridSpec: class header, task header, heatmap, db labels, colorbar
    gs = GridSpec(5, n_width_cols, figure=fig,
                  height_ratios=[0.4, task_header_height, 10, 1.5, 0.3],
                  width_ratios=width_ratios,
                  wspace=0.02, hspace=0.02)

    # --- Method names (leftmost column, spanning heatmap rows) ---
    ax_names = fig.add_subplot(gs[2, 0])
    ax_names.set_xlim(0, 1)
    ax_names.set_ylim(-0.5, n_methods - 0.5)
    ax_names.invert_yaxis()
    ax_names.axis('off')

    for i, method in enumerate(method_order):
        fontweight = 'normal' if method in baselines else 'bold'
        ax_names.text(1, i, method, ha='right', va='center',
                     fontsize=9, color='black', fontweight=fontweight)

    # Shared colormap settings
    vmin = 0
    vmax = np.nanmax(db_mean.values)
    cmap = plt.cm.viridis
    threshold = (vmin + vmax) / 2

    im = None

    # Track class spans for class headers
    class_spans = {}  # {cls: (start_gs_col, end_gs_col)}
    for i, (cls, task, cols) in enumerate(task_list):
        gs_col = gs_col_mapping[i]
        if cls not in class_spans:
            class_spans[cls] = [gs_col, gs_col]
        else:
            class_spans[cls][1] = gs_col

    # --- Class headers (spanning multiple task columns) ---
    for cls, (start_col, end_col) in class_spans.items():
        ax_class = fig.add_subplot(gs[0, start_col:end_col + 1])
        ax_class.set_xlim(0, 1)
        ax_class.set_ylim(0, 1)
        ax_class.axis('off')
        rect = plt.Rectangle((0, 0), 1, 1,
                              facecolor=class_colors.get(cls, 'white'),
                              edgecolor='black', linewidth=1)
        ax_class.add_patch(rect)
        label = class_names_map.get(cls, cls)
        ax_class.text(0.5, 0.5, label, ha='center', va='center',
                      fontsize=10, fontweight='bold')

    # --- Render each task as separate heatmap ---
    for i, (cls, task, cols) in enumerate(task_list):
        gs_col = gs_col_mapping[i]
        n_task_cols = len(cols)

        # --- Task header ---
        ax_task = fig.add_subplot(gs[1, gs_col])
        ax_task.set_xlim(0, 1)
        ax_task.set_ylim(0, 1)
        ax_task.axis('off')
        rect = plt.Rectangle((0, 0), 1, 1,
                              facecolor='white', edgecolor='gray', linewidth=0.5)
        ax_task.add_patch(rect)
        wrapped_label = task_wrapped_titles[task]
        ax_task.text(0.5, 0.5, wrapped_label, ha='center', va='center', fontsize=9)

        # --- Heatmap for this task ---
        ax_heatmap = fig.add_subplot(gs[2, gs_col])
        heatmap_data = db_mean[cols].values

        im = ax_heatmap.imshow(heatmap_data, aspect='auto', cmap=cmap, vmin=vmin, vmax=vmax)

        # Add ranking text
        for row_i in range(n_methods):
            for col_j in range(n_task_cols):
                col = cols[col_j]
                rank_val = db_ranks.loc[method_order[row_i], col]
                cell_val = heatmap_data[row_i, col_j]
                if pd.isna(rank_val) or pd.isna(cell_val):
                    rank_str = '-'
                    text_color = 'gray'
                else:
                    text_color = 'white' if cell_val > threshold else 'black'
                    if rank_val == int(rank_val):
                        rank_str = f'{int(rank_val)}'
                    else:
                        rank_str = f'{rank_val:.1f}'
                ax_heatmap.text(col_j, row_i, rank_str, ha='center', va='center',
                               fontsize=10, color=text_color)

        ax_heatmap.set_yticks([])
        ax_heatmap.set_xticks([])

        # --- Database labels at bottom (rotated 90 degrees) ---
        ax_db = fig.add_subplot(gs[3, gs_col])
        ax_db.set_xlim(-0.5, n_task_cols - 0.5)
        ax_db.set_ylim(0, 1)
        ax_db.axis('off')

        for j, (db, _) in enumerate(cols):
            label = db_names_map.get(db, db)
            ax_db.text(j, 1, label, ha='center', va='top', fontsize=10, rotation=90)

    # --- Colorbar (small, centered at bottom) ---
    cbar_width = 0.15
    cbar_left = 0.5 - cbar_width / 2
    cbar_ax = fig.add_axes([cbar_left, 0.02, cbar_width, 0.015])
    cbar = plt.colorbar(im, cax=cbar_ax, orientation='horizontal')
    cbar.set_label('Mean F0.1', fontsize=9)
    cbar.ax.tick_params(labelsize=7)

    return fig


def create_pair_comparison_figure(pair_df, config):
    """Create barplot figure comparing Unpaired vs Paired and Synthetic vs Paired."""
    # Get method colors from config with name mapping
    method_colors_old = config['colors']['nets']
    method_colors = {config['method_names'][k]: method_colors_old[k] for k in method_colors_old}

    # Filter data into two subsets (note: typo "Synthethic" in CSV)
    unpaired_df = pair_df[pair_df['type'] == 'Unpaired vs Paired'].copy()
    synthetic_df = pair_df[pair_df['type'] == 'Synthethic vs Paired'].copy()

    # Order methods by f01 values from "Unpaired vs Paired" comparison
    unpaired_df = unpaired_df.sort_values('f01', ascending=False)
    method_order = unpaired_df['name'].tolist()

    # Reorder synthetic_df to match
    synthetic_df = synthetic_df.set_index('name').loc[method_order].reset_index()

    # Create 1x2 subplot figure
    fig, axes = plt.subplots(1, 2, figsize=(7, 3.5), sharex=True, sharey=True)

    # Plot Unpaired vs Paired (left panel)
    ax1 = axes[0]
    bar_colors = [method_colors.get(m, 'gray') for m in method_order]
    y_positions = np.arange(len(method_order))
    ax1.barh(y_positions, unpaired_df['f01'].values, color=bar_colors, height=0.8)
    ax1.axvline(x=0, ls='--', lw=1, c='black')
    ax1.set_yticks(y_positions)
    ax1.set_yticklabels(method_order, fontsize=9)
    ax1.set_xlabel('ΔF0.1', fontsize=10)
    ax1.set_title('Unpaired vs Paired', fontsize=11)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # Plot Synthetic vs Paired (right panel)
    ax2 = axes[1]
    ax2.barh(y_positions, synthetic_df['f01'].values, color=bar_colors, height=0.8)
    ax2.axvline(x=0, ls='--', lw=1, c='black')
    ax2.set_yticks(y_positions)
    ax2.set_yticklabels(method_order, fontsize=9)
    ax2.set_xlabel('ΔF0.1', fontsize=10)
    ax2.set_title('Synthetic vs Paired', fontsize=11)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    plt.tight_layout()
    return fig


def create_topology_correlation_figure(df, stats_df, config):
    """Create scatter plots of network topology statistics vs F0.1 performance.

    Args:
        df: Metrics DataFrame with columns 'dts', 'name', 'f01'
        stats_df: Topology stats DataFrame with columns 'dts', 'name', '# TFs', etc.
        config: Configuration dictionary with method colors

    Returns:
        matplotlib Figure with 5 scatter plots and shared legend
    """
    # Get method colors from config
    method_colors_old = config['colors']['nets']
    method_colors = {config['method_names'][k]: method_colors_old[k] for k in method_colors_old}

    # Merge metrics and stats on ['dts', 'name']
    merged = df.merge(stats_df, on=['dts', 'name'], how='inner')

    # Aggregate by method name: compute mean of numeric columns
    agg_df = merged.groupby('name').mean(numeric_only=True).reset_index()

    # Stat columns to plot
    stat_columns = ['# TFs', '# CREs', '# Genes', '# Edges', 'Regulon size']

    # Create figure: 1 row x 5 columns
    fig, axes = plt.subplots(1, 5, figsize=(12, 3), sharey=True)

    # Get all method names for consistent legend
    method_names = agg_df['name'].unique()

    # First pass: compute all correlations and p-values
    correlations = []
    for i, stat_col in enumerate(stat_columns):
        x = agg_df[stat_col].values
        y = agg_df['f01'].values
        mask = ~(np.isnan(x) | np.isnan(y))
        if mask.sum() > 2:
            r, p = pearsonr(x[mask], y[mask])
            correlations.append({'stat': stat_col, 'r': r, 'p': p})
        else:
            correlations.append({'stat': stat_col, 'r': np.nan, 'p': np.nan})

    # Apply FDR correction (Benjamini-Hochberg)
    p_values = np.array([c['p'] for c in correlations])
    valid_mask = ~np.isnan(p_values)
    if valid_mask.sum() > 0:
        fdr_adjusted = benjamini_hochberg(p_values[valid_mask])
        fdr_full = np.full(len(p_values), np.nan)
        fdr_full[valid_mask] = fdr_adjusted
        for i, c in enumerate(correlations):
            c['fdr'] = fdr_full[i]

    # Second pass: create plots with FDR-corrected significance
    for i, stat_col in enumerate(stat_columns):
        ax = axes[i]

        # Get data for this stat
        x = agg_df[stat_col].values
        y = agg_df['f01'].values
        names = agg_df['name'].values

        # Scatter plot colored by method
        for j, name in enumerate(names):
            color = method_colors.get(name, 'gray')
            ax.scatter(x[j], y[j], c=color, s=50, label=name if i == 0 else None,
                      edgecolors='white', linewidth=0.5)

        # Get correlation result for this stat
        corr = correlations[i]
        if not np.isnan(corr['r']):
            r = corr['r']
            fdr = corr['fdr']

            # Significance stars based on FDR < 0.05
            if fdr < 0.001:
                stars = '***'
            elif fdr < 0.01:
                stars = '**'
            elif fdr < 0.05:
                stars = '*'
            else:
                stars = ''

            # Annotate with r value and stars (bottom right)
            ax.text(0.95, 0.05, f'r = {r:.2f} {stars}', transform=ax.transAxes,
                   fontsize=9, va='bottom', ha='right')

        ax.set_xlabel(stat_col, fontsize=9)
        if i == 0:
            ax.set_ylabel('Mean F0.1', fontsize=9)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(labelsize=8)

    # Create shared legend below the plots
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 0),
              ncol=5, fontsize=8, frameon=False)

    # Adjust layout to make room for legend
    fig.subplots_adjust(bottom=0.35)

    return fig


def main():
    """Main function to generate the ranking figure."""
    # Parse arguments
    parser = argparse.ArgumentParser(description='Generate GRN benchmark ranking figure')
    parser.add_argument('-i', '--path_inp', required=True)
    parser.add_argument('-s', '--path_scalability', required=True,
                        help='Path to scalability CSV (columns: mth, h, gb, use_gpu)')
    parser.add_argument('-p', '--path_pair', required=True,
                        help='Path to pair comparison CSV')
    parser.add_argument('-t', '--path_stats', required=True,
                        help='Path to topology stats CSV (columns: dts, name, # TFs, etc.)')
    parser.add_argument('-o', '--path_out', required=True)
    args = vars(parser.parse_args())
    path_inp = args['path_inp']
    path_scalability = args['path_scalability']
    path_pair = args['path_pair']
    path_stats = args['path_stats']
    path_out = args['path_out']

    # Load data
    df, config = load_data(path_inp)

    # Load pair data
    pair_df = load_pair_data(path_pair)

    # Load scalability data
    scalability_df = load_scalability_data(path_scalability)

    # Load topology stats data
    stats_df = load_stats_data(path_stats)

    # Compute aggregations
    overall_mean, class_mean, dataset_mean = compute_aggregations(df)

    # Compute rankings
    class_ranks = compute_rankings(class_mean)
    dataset_ranks = compute_rankings(dataset_mean)

    # Compute scalability rankings (lower is better)
    scalability_ranks = compute_scalability_rankings(scalability_df)

    # Compute mean rank across the 4 classes
    mean_class_rank = class_ranks.mean(axis=1)

    # Sort datasets by mean performance across methods (descending) for consistent ordering
    dataset_order = dataset_mean.mean().sort_values(ascending=False).index.tolist()

    # Compute dataset comparison statistics
    kw_stat, kw_pval, posthoc_df = compute_dataset_comparison_stats(df)

    # Create ranking figure (page 1)
    fig1 = create_figure(overall_mean, class_mean, dataset_mean,
                         class_ranks, dataset_ranks,
                         scalability_df, scalability_ranks, config,
                         use_rank_barplot=True,
                         mean_class_rank=mean_class_rank)

    # Create dataset comparison figure (page 2)
    fig2 = create_dataset_boxplot_figure(dataset_mean, posthoc_df, config, dataset_order)

    # Compute database-level aggregation (per db-task combination)
    db_mean, db_hierarchy = compute_database_aggregation(df)

    # Order (db, task) combinations hierarchically
    # Use full class names from data, ordered by config's class_names keys
    class_names = config.get('class_names', {})
    db_class_order = [class_names.get(c, c) for c in ['genom', 'pred', 'prior', 'mech']]
    ordered_cols = get_hierarchical_db_order(db_mean, db_hierarchy, db_class_order)

    # Reorder and compute rankings
    db_mean = db_mean[ordered_cols]
    db_ranks = compute_rankings(db_mean)

    # Get method order (same as main figure - sorted by mean rank ascending)
    method_order = mean_class_rank.sort_values(ascending=True).index.tolist()

    # Create database heatmap figure (page 3)
    fig3 = create_database_heatmap_figure(
        db_mean.loc[method_order],
        db_ranks.loc[method_order],
        db_hierarchy,
        method_order,
        db_class_order,
        config
    )

    # Create pair comparison figure (page 4)
    fig4 = create_pair_comparison_figure(pair_df, config)

    # Create topology correlation figure (page 5)
    fig5 = create_topology_correlation_figure(df, stats_df, config)

    # Save multi-page PDF
    with PdfPages(path_out) as pdf:
        pdf.savefig(fig1, bbox_inches='tight', dpi=300)
        pdf.savefig(fig2, bbox_inches='tight', dpi=300)
        pdf.savefig(fig3, bbox_inches='tight', dpi=300)
        pdf.savefig(fig4, bbox_inches='tight', dpi=300)
        pdf.savefig(fig5, bbox_inches='tight', dpi=300)

    plt.close('all')


if __name__ == '__main__':
    main()
