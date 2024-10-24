import sys
import os
import warnings
os.environ['NUMBA_CACHE_DIR'] = '/tmp/'
import scanpy as sc
import anndata as ad
from typing import (List, Optional, Tuple, Set, Dict, Union, Literal, TYPE_CHECKING)
import pathlib
from pathlib import Path
import argparse
import logging
import pickle
import requests

import joblib
import numpy as np
import pandas as pd
import h5py
import scipy as sp
from scipy.stats import gaussian_kde
import muon as mu

from matplotlib import cm
import subprocess
from matplotlib.colors import Normalize
from scipy.stats import pearsonr, spearmanr
from sklearn.ensemble import (ExtraTreesRegressor, GradientBoostingRegressor,
                              RandomForestRegressor)
from tqdm import tqdm

# for chromsizes
import polars as pl
pl.enable_string_cache()

import pyranges as pr
import pysam
import pycisTopic.pseudobulk_peak_calling
from pycisTopic.fragments import (
    get_fragments_in_peaks,
    get_fragments_per_cb,
    get_insert_size_distribution,
)
from pycisTopic.topic_binarization import threshold_otsu, binarize_topics
from pycisTopic.tss_profile import get_tss_profile
from pycisTopic.iterative_peak_calling import get_consensus_peaks
from pycisTopic.qc import get_barcodes_passing_qc_for_sample
from pycisTopic.cistopic_class import create_cistopic_object_from_fragments
import pycisTopic.lda_models
import pycisTopic.cistopic_class
from pycisTopic.diff_features import (
    impute_accessibility,
    normalize_scores,
    find_highly_variable_features,
    find_diff_features)
from pycisTopic.utils import region_names_to_coordinates

import scenicplus.data_wrangling.adata_cistopic_wrangling
import scenicplus.data_wrangling.gene_search_space

import pycistarget.motif_enrichment_cistarget
from pycistarget.motif_enrichment_dem import (
    DEM,
    ranksums_numba_multiple,
    mean_axis1, get_log2_fc,
    p_adjust_bh, get_optimal_threshold_roc,
    DEMDatabase
    )
from scenicplus.triplet_score import get_max_rank_of_motif_for_each_TF
import mudata
from scenicplus.cli.commands import (
#    prepare_motif_enrichment_results,
    _get_foreground_background
)

from pycistarget.motif_enrichment_result import MotifEnrichmentResult

from scenicplus.utils import p_adjust_bh
from scenicplus.grn_builder.modules import (
    create_emodules, eRegulon, merge_emodules, RHO_THRESHOLD, TARGET_GENE_NAME)
from scenicplus.grn_builder.gsea import run_gsea
from scenicplus.grn_builder.gsea_approach import _run_gsea_for_e_module
from scenicplus.triplet_score import (
        get_max_rank_of_motif_for_each_TF,
        _rank_scores_and_assign_random_ranking_in_range_for_ties,
        _calculate_cross_species_rank_ratio_with_order_statistics)
import anndata
from tqdm import tqdm


log = logging.getLogger("SCENIC+")

score_to_keep = "rho"  #"importance_x_rho"

RANDOM_SEED = 666

SKLEARN_REGRESSOR_FACTORY = {
    'RF': RandomForestRegressor,
    'ET': ExtraTreesRegressor,
    'GBM': GradientBoostingRegressor
}

SCIPY_CORRELATION_FACTORY = {
    'PR': pearsonr,
    'SR': spearmanr
}

# Parameters from arboreto
# scikit-learn random forest regressor
RF_KWARGS = {
    'n_jobs': 1,
    'n_estimators': 1000,
    'max_features': 'sqrt'
}

# scikit-learn extra-trees regressor
ET_KWARGS = {
    'n_jobs': 1,
    'n_estimators': 1000,
    'max_features': 'sqrt'
}

# scikit-learn gradient boosting regressor
GBM_KWARGS = {
    'learning_rate': 0.01,
    'n_estimators': 500,
    'max_features': 0.1
}

# scikit-learn stochastic gradient boosting regressor
SGBM_KWARGS = {
    'learning_rate': 0.01,
    'n_estimators': 5000,  # can be arbitrarily large
    'max_features': 0.1,
    'subsample': 0.9
}

def compute_qc_stats(
    fragments_df_pl: pl.DataFrame,
    regions_df_pl: pl.DataFrame,
    tss_annotation: pl.DataFrame,
    tss_flank_window: int = 2000,
    tss_smoothing_rolling_window: int = 10,
    tss_minimum_signal_window: int = 100,
    tss_window: int = 50,
    tss_min_norm: float = 0.2,
    use_genomic_ranges: bool = True,
    min_fragments_per_cb: int = 10,
    collapse_duplicates: bool = True,
    no_threads: int = 8,
) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    """
    Compute quality check statistics from Polars DataFrame with fragments.

    Parameters
    ----------
    fragments_df_pl
        Polars DataFrame with fragments.
        fragments_df_pl
        Polars DataFrame with fragments (filtered by cell barcodes of interest).
        See :func:`pycisTopic.fragments.filter_fragments_by_cb`.
    regions_df_pl
        Polars DataFrame with peak regions (consensus peaks or SCREEN regions).
        See :func:`pycisTopic.fragments.read_bed_to_polars_df` for a way to read a BED
        file with peak regions.
    tss_annotation
        TSS annotation Polars DataFrame with at least the following columns:
        ``["Chromosome", "Start", "Strand"]``.
        The "Start" column is 0-based like a BED file.
        See :func:`pycisTopic.gene_annotation.read_tss_annotation_from_bed`,
        :func:`pycisTopic.gene_annotation.get_tss_annotation_from_ensembl` and
        :func:`pycisTopic.gene_annotation.change_chromosome_source_in_bed` for ways
        to get TSS annotation from Ensembl BioMart.
    tss_flank_window
        Flanking window around the TSS.
        Used for intersecting fragments with TSS positions and keeping cut sites.
        Default: ``2000`` (+/- 2000 bp).
        See :func:`pycisTopic.tss_profile.get_tss_profile`.
    tss_smoothing_rolling_window
        Rolling window used to smooth the cut sites signal.
        Default: ``10``.
        See :func:`pycisTopic.tss_profile.get_tss_profile`.
    tss_minimum_signal_window
        Average signal in the tails of the flanking window around the TSS:
           - ``[-flank_window, -flank_window + minimum_signal_window + 1]``
           - ``[flank_window - minimum_signal_window + 1, flank_window]``
        is used to normalize the TSS enrichment.
        Default: ``100`` (average signal in ``[-2000, -1901]``, ``[1901, 2000]``
        around TSS if ``flank_window=2000``).
        See :func:`pycisTopic.tss_profile.get_tss_profile`.
    tss_window
        Window around the TSS used to count fragments in the TSS when calculating
        the TSS enrichment per cell barcode.
        Default: ``50`` (+/- 50 bp).
        See :func:`pycisTopic.tss_profile.get_tss_profile`.
    tss_min_norm
        Minimum normalization score.
        If the average minimum signal value is below this value, this number is used
        to normalize the TSS signal. This approach penalizes cells with fewer reads.
        Default: ``0.2``
        See :func:`pycisTopic.tss_profile.get_tss_profile`.
    use_genomic_ranges
        Use genomic ranges implementation for calculating intersections, instead of
        using pyranges.
    min_fragments_per_cb
        Minimum number of fragments needed per cell barcode to keep the fragments
        for those cell barcodes.
    collapse_duplicates
        Collapse duplicate fragments (same chromosomal positions and linked to the same
        cell barcode).
    no_threads
        Number of threads to use when calculating kernel-density estimate (KDE) to get
        probability density function (PDF) values for log10 unique fragments in peaks
        vs TSS enrichment, fractions of fragments in peaks and duplication ratio.
        Default: ``8``

    Returns
    -------
    Tuple with:
      - Polars DataFrame with fragments statistics per cell barcode.
      - Polars DataFrame with insert size distribution of fragments.
      - Polars DataFrame with TSS normalization matrix for the whole sample.
      - Polars DataFrame with TSS normalization matrix per cell barcode.

    See Also
    --------
    pycisTopic.fragments.filter_fragments_by_cb
    pycisTopic.fragments.get_insert_size_distribution
    pycisTopic.fragments.get_fragments_in_peaks
    pycisTopic.fragments.read_bed_to_polars_df
    pycisTopic.fragments.read_fragments_to_polars_df
    pycisTopic.gene_annotation.read_tss_annotation_from_bed
    pycisTopic.tss_profile.get_tss_profile

    Examples
    --------
    >>> from pycisTopic.fragments import read_bed_to_polars_df
    >>> from pycisTopic.fragments import read_fragments_to_polars_df
    >>> from pycisTopic.gene_annotation import read_tss_annotation_from_bed

    1. Read gzipped fragments BED file to a Polars DataFrame.

    >>> fragments_df_pl = read_fragments_to_polars_df(
    ...     fragments_bed_filename="fragments.tsv.gz",
    ... )

    2. Read BED file with consensus peaks or SCREEN regions (get first 3 columns only)
       which will be used for counting number of fragments in peaks.

    >>> regions_df_pl = read_bed_to_polars_df(
    ...     bed_filename=screen_regions_bed_filename,
    ...     min_column_count=3,
    ... )

    3. Read TSS annotation from a file.
       See :func:`pycisTopic.gene_annotation.read_tss_annotation_from_bed` for more
       info.

    >>> tss_annotation_bed_df_pl = read_tss_annotation_from_bed(
    ...     tss_annotation_bed_filename="hg38.tss.bed",
    ... )

    4. Compute QC statistics.

    >>> (
    ...     fragments_stats_per_cb_df_pl,
    ...     insert_size_dist_df_pl,
    ...     tss_norm_matrix_sample,
    ...     tss_norm_matrix_per_cb,
    ... ) = compute_qc_stats(
    ...     fragments_df_pl=fragments_cb_filtered_df_pl,
    ...     regions_df_pl=regions_df_pl,
    ...     tss_annotation=tss_annotation_bed_df_pl,
    ...     tss_flank_window=2000,
    ...     tss_smoothing_rolling_window=10,
    ...     tss_minimum_signal_window=100,
    ...     tss_window=50,
    ...     tss_min_norm=0.2,
    ...     use_genomic_ranges=True,
    ...     min_fragments_per_cb=10,
    ...     collapse_duplicates=True,
    ...     no_threads=8,
    ... )

    """
    from pycisTopic.qc import compute_kde
    logger = logging.getLogger(__name__)
    # Define correct column to get, based on the setting of `collapse_duplicates`.
    fragments_count_column = (
        "unique_fragments_count" if collapse_duplicates else "total_fragments_count"
    )
    fragments_in_peaks_count_column = (
        "unique_fragments_in_peaks_count"
        if collapse_duplicates
        else "total_fragments_in_peaks_count"
    )
    # Get Polars DataFrame with basic fragments statistics per cell barcode.
    logger.info("Get basic fragments statistics per cell barcode.")
    fragments_stats_per_cb_df_pl = get_fragments_per_cb(
        fragments_df_pl=fragments_df_pl,
        min_fragments_per_cb=min_fragments_per_cb,
        collapse_duplicates=collapse_duplicates,
    ).rename({"by":"CB"})
    # Get Polars DataFrame with total fragment counts and unique fragment counts
    # per region.
    logger.info("Get total fragment counts and unique fragment counts per region.")
    fragments_in_peaks_df_pl = get_fragments_in_peaks(
        fragments_df_pl=fragments_df_pl,
        regions_df_pl=regions_df_pl,
    ).rename({"by": "CB"})
    # Add fragment counts per region to fragments statistics per cell barcode.
    logger.info(
        "Add fragment counts per region to fragments statistics per cell barcode."
    )
    fragments_stats_per_cb_df_pl = (
        fragments_stats_per_cb_df_pl.lazy()
        .join(
            fragments_in_peaks_df_pl.lazy(),
            how="left",
            on="CB",
        )
        .with_columns(
            pl.col("total_fragments_in_peaks_count").fill_null(0),
            pl.col("unique_fragments_in_peaks_count").fill_null(0),
        )
        .with_columns(
            (
                pl.col(fragments_in_peaks_count_column) / pl.col(fragments_count_column)
            ).alias("fraction_of_fragments_in_peaks")
        )
        .select(
            pl.col("CB"),
            pl.col("barcode_rank"),
            pl.col("total_fragments_count"),
            (pl.col("total_fragments_count") + 1)
            .log10()
            .alias("log10_total_fragments_count"),
            pl.col("unique_fragments_count"),
            (pl.col("unique_fragments_count") + 1)
            .log10()
            .alias("log10_unique_fragments_count"),
            pl.col("total_fragments_in_peaks_count"),
            (pl.col("total_fragments_in_peaks_count") + 1)
            .log10()
            .alias("log10_total_fragments_in_peaks_count"),
            pl.col("unique_fragments_in_peaks_count"),
            (pl.col("unique_fragments_in_peaks_count") + 1)
            .log10()
            .alias("log10_unique_fragments_in_peaks_count"),
            pl.col("fraction_of_fragments_in_peaks"),
            pl.col("duplication_count"),
            pl.col("duplication_ratio"),
            pl.col("nucleosome_signal"),
        )
    )
    # Get insert size distribution of fragments.
    logger.info("Get insert size distribution of fragments.")
    insert_size_dist_df_pl = get_insert_size_distribution(
        fragments_df_pl=fragments_df_pl,
    )
    # Get TSS profile for fragments.
    logger.info("Get TSS profile for fragments.")
    (
        tss_enrichment_per_cb,
        tss_norm_matrix_sample,
        tss_norm_matrix_per_cb,
    ) = get_tss_profile(
        fragments_df_pl=fragments_df_pl,
        tss_annotation=tss_annotation,
        flank_window=tss_flank_window,
        smoothing_rolling_window=tss_smoothing_rolling_window,
        minimum_signal_window=tss_minimum_signal_window,
        tss_window=tss_window,
        min_norm=tss_min_norm,
        use_genomic_ranges=use_genomic_ranges,
    )
    # Add TSS enrichment to fragments statistics per cell barcode.
    logger.info("Add TSS enrichment to fragments statistics per cell barcode.")
    print(fragments_stats_per_cb_df_pl)
    print(tss_enrichment_per_cb)
    fragments_stats_per_cb_df_pl = (
        fragments_stats_per_cb_df_pl.join(
            tss_enrichment_per_cb.lazy(),
            how="left",
            on="CB",
        )
        .with_columns(
            pl.col("tss_enrichment").fill_null(0.0),
        )
        .collect()
    )
    # Extract certain columns as numpy arrays as they are needed for calculating KDE.
    (
        log10_unique_fragments_in_peaks_count,
        tss_enrichment,
        fraction_of_fragments_in_peaks,
        duplication_ratio,
    ) = (
        fragments_stats_per_cb_df_pl.select(
            [
                pl.col("log10_unique_fragments_in_peaks_count"),
                pl.col("tss_enrichment"),
                pl.col("fraction_of_fragments_in_peaks"),
                pl.col("duplication_ratio"),
            ]
        )
        .to_numpy()
        .T
    )
    # Construct 2D numpy matrices for usage with compute_kde.
    kde_data_for_tss_enrichment = np.vstack(
        [log10_unique_fragments_in_peaks_count, tss_enrichment]
    )
    kde_data_for_fraction_of_fragments_in_peaks = np.vstack(
        [log10_unique_fragments_in_peaks_count, fraction_of_fragments_in_peaks]
    )
    kde_data_for_duplication_ratio = np.vstack(
        [log10_unique_fragments_in_peaks_count, duplication_ratio]
    )
    # Calculate KDE for log10 unique fragments in peaks vs TSS enrichment,
    # fractions of fragments in peaks and duplication ratio.
    logger.info("Calculate KDE for log10 unique fragments in peaks vs TSS enrichment.")
    pdf_values_for_tss_enrichment = compute_kde(
        training_data=kde_data_for_tss_enrichment,
        test_data=kde_data_for_tss_enrichment,
        no_threads=no_threads,
    )
    logger.info(
        "Calculate KDE for log10 unique fragments in peaks vs fractions of fragments "
        "in peaks."
    )
    pdf_values_for_fraction_of_fragments_in_peaks = compute_kde(
        training_data=kde_data_for_fraction_of_fragments_in_peaks,
        test_data=kde_data_for_fraction_of_fragments_in_peaks,
        no_threads=no_threads,
    )
    logger.info(
        "Calculate KDE for log10 unique fragments in peaks vs duplication ratio."
    )
    pdf_values_for_duplication_ratio = compute_kde(
        training_data=kde_data_for_duplication_ratio,
        test_data=kde_data_for_duplication_ratio,
        no_threads=no_threads,
    )
    # Add probability density function (PDF) values for log10 unique fragments in peaks
    # vs TSS enrichment, fractions of fragments in peaks and duplication ratio to
    # fragments statistics per cell barcode.
    logger.info(
        "Add probability density function (PDF) values to fragments statistics per "
        "cell barcode."
    )
    fragments_stats_per_cb_df_pl = fragments_stats_per_cb_df_pl.hstack(
        pl.DataFrame(
            {
                "pdf_values_for_tss_enrichment": pdf_values_for_tss_enrichment,
                "pdf_values_for_fraction_of_fragments_in_peaks": pdf_values_for_fraction_of_fragments_in_peaks,
                "pdf_values_for_duplication_ratio": pdf_values_for_duplication_ratio,
            }
        )
    )
    return (
        fragments_stats_per_cb_df_pl,
        insert_size_dist_df_pl,
        tss_norm_matrix_sample,
        tss_norm_matrix_per_cb,
    )


def qc(
    fragments_tsv_filename: str | Path,
    regions_bed_filename: str | Path | pl.DataFrame,
    tss_annotation_bed_df_pl: pl.DataFrame,
    output_prefix: str,
    tss_flank_window: int = 2000,
    tss_smoothing_rolling_window: int = 10,
    tss_minimum_signal_window: int = 100,
    tss_window: int = 50,
    tss_min_norm: float = 0.2,
    use_genomic_ranges: bool = True,
    min_fragments_per_cb: int = 10,
    collapse_duplicates: bool = True,
    no_threads: int = 8,
    engine: str | Literal["polars"] | Literal["pyarrow"] = "pyarrow",
) -> None:
    """
    Compute quality check statistics from fragments file.

    Parameters
    ----------
    fragments_tsv_filename
        Fragments TSV filename which contains scATAC fragments.
    regions_bed_filename
        Consensus peaks / SCREEN regions BED file.
        Used to calculate amount of fragments in peaks.
        Can be a string/Path to a BED file or a Polars DataFrame. 
    tss_annotation_bed_filename
        TSS annotation BED file.
        Used to calculate distance of fragments to TSS positions.
    output_prefix
        Output prefix to use for QC statistics parquet output files.
    tss_flank_window
        Flanking window around the TSS.
        Used for intersecting fragments with TSS positions and keeping cut sites.
        Default: ``2000`` (+/- 2000 bp).
        See :func:`pycisTopic.tss_profile.get_tss_profile`.
    tss_smoothing_rolling_window
        Rolling window used to smooth the cut sites signal.
        Default: ``10``.
        See :func:`pycisTopic.tss_profile.get_tss_profile`.
    tss_minimum_signal_window
        Average signal in the tails of the flanking window around the TSS:
           - ``[-flank_window, -flank_window + minimum_signal_window + 1]``
           - ``[flank_window - minimum_signal_window + 1, flank_window]``
        is used to normalize the TSS enrichment.
        Default: ``100`` (average signal in ``[-2000, -1901]``, ``[1901, 2000]``
        around TSS if `flank_window=2000`).
        See :func:`pycisTopic.tss_profile.get_tss_profile`.
    tss_window
        Window around the TSS used to count fragments in the TSS when calculating
        the TSS enrichment per cell barcode.
        Default: ``50`` (+/- 50 bp).
        See :func:`pycisTopic.tss_profile.get_tss_profile`.
    tss_min_norm
        Minimum normalization score.
        If the average minimum signal value is below this value, this number is used
        to normalize the TSS signal. This approach penalizes cells with fewer reads.
        Default: ``0.2``
        See :func:`pycisTopic.tss_profile.get_tss_profile`.
    use_genomic_ranges
        Use genomic ranges implementation for calculating intersections, instead of
        using pyranges.
    min_fragments_per_cb
        Minimum number of fragments needed per cell barcode to keep the fragments
        for those cell barcodes.
    collapse_duplicates
        Collapse duplicate fragments (same chromosomal positions and linked to the same
        cell barcode).
    no_threads
        Number of threads to use when calculating kernel-density estimate (KDE) to get
        probability density function (PDF) values for log10 unique fragments in peaks
        vs TSS enrichment, fractions of fragments in peaks and duplication ratio.
        Default: ``8``
    engine
        Use Polars or pyarrow to read BED and fragment files (default: `pyarrow`).

    Returns
    -------
    None

    """
    from pycisTopic.fragments import read_bed_to_polars_df, read_fragments_to_polars_df
    from pycisTopic.gene_annotation import read_tss_annotation_from_bed
    from pycisTopic.qc import get_otsu_threshold
    # Remove trailing dot(s) from the output prefix.
    output_prefix = output_prefix.rstrip(".")
    if isinstance(regions_bed_filename, pl.DataFrame):
        regions_df_pl = regions_bed_filename
    else:
        print(f'Loading regions BED file from "{regions_bed_filename}".')
        regions_df_pl = read_bed_to_polars_df(
            bed_filename=regions_bed_filename,
            min_column_count=3,
            engine=engine,
        )
    print(f'Loading fragments TSV file from "{fragments_tsv_filename}".')
    fragments_df_pl = read_fragments_to_polars_df(
        fragments_tsv_filename,
        engine=engine,
    )
    print("Computing QC stats.")
    (
        fragments_stats_per_cb_df_pl,
        insert_size_dist_df_pl,
        tss_norm_matrix_sample,
        tss_norm_matrix_per_cb,
    ) = compute_qc_stats(
        fragments_df_pl=fragments_df_pl,
        regions_df_pl=regions_df_pl,
        tss_annotation=tss_annotation_bed_df_pl,
        tss_flank_window=tss_flank_window,
        tss_smoothing_rolling_window=tss_smoothing_rolling_window,
        tss_minimum_signal_window=tss_minimum_signal_window,
        tss_window=tss_window,
        tss_min_norm=tss_min_norm,
        use_genomic_ranges=use_genomic_ranges,
        min_fragments_per_cb=min_fragments_per_cb,
        collapse_duplicates=collapse_duplicates,
        no_threads=no_threads,
    )
    print(f'Writing "{output_prefix}.fragments_stats_per_cb.parquet".')
    fragments_stats_per_cb_df_pl.write_parquet(
        f"{output_prefix}.fragments_stats_per_cb.parquet",
        compression="zstd",
        use_pyarrow=True,
    )
    print(f'Writing "{output_prefix}.fragments_insert_size_dist.parquet".')
    insert_size_dist_df_pl.write_parquet(
        f"{output_prefix}.fragments_insert_size_dist.parquet",
        compression="zstd",
        use_pyarrow=True,
    )
    print(f'Writing "{output_prefix}.tss_norm_matrix_sample.parquet".')
    tss_norm_matrix_sample.write_parquet(
        f"{output_prefix}.tss_norm_matrix_sample.parquet",
        compression="zstd",
        use_pyarrow=True,
    )
    print(f'Writing "{output_prefix}.tss_norm_matrix_per_cb.parquet".')
    tss_norm_matrix_per_cb.write_parquet(
        f"{output_prefix}.tss_norm_matrix_per_cb.parquet",
        compression="zstd",
        use_pyarrow=True,
    )
    print("Calculating Otsu thresholds.")
    (
        unique_fragments_in_peaks_count_otsu_threshold,
        tss_enrichment_otsu_threshold,
        fragments_stats_per_cb_for_otsu_threshold_df_pl,
    ) = get_otsu_threshold(
        fragments_stats_per_cb_df_pl=fragments_stats_per_cb_df_pl,
        min_otsu_fragments=100,
        min_otsu_tss=1.0,
    )
    print(
        f'Writing "{output_prefix}.fragments_stats_per_cb_for_otsu_thresholds.parquet".'
    )
    fragments_stats_per_cb_for_otsu_threshold_df_pl.write_parquet(
        f"{output_prefix}.fragments_stats_per_cb_for_otsu_thresholds.parquet",
        compression="zstd",
        use_pyarrow=True,
    )
    print(
        f'Writing "{output_prefix}.fragments_stats_per_cb_for_otsu_thresholds.tsv".'
    )
    fragments_stats_per_cb_for_otsu_threshold_df_pl.write_csv(
        f"{output_prefix}.fragments_stats_per_cb_for_otsu_thresholds.tsv",
        separator="\t",
        include_header=True,
    )
    print(f'Writing "{output_prefix}.cbs_for_otsu_thresholds.tsv".')
    fragments_stats_per_cb_for_otsu_threshold_df_pl.select(
        pl.col("CB")).write_csv(
            f"{output_prefix}.cbs_for_otsu_thresholds.tsv",
            separator="\t",
            include_header=False,
            )
    print(f'Writing "{output_prefix}.otsu_thresholds.tsv".')
    with open(f"{output_prefix}.otsu_thresholds.tsv", "w") as fh:
        print(
            "unique_fragments_in_peaks_count_otsu_threshold\t" +
            "tss_enrichment_otsu_threshold\n"
            f"{unique_fragments_in_peaks_count_otsu_threshold}\t{tss_enrichment_otsu_threshold}",
            file=fh,
        )
    print("pycisTopic QC finished.")


def create_fragment_matrix_from_fragments(
    fragments_bed_filename: str | Path,
    regions_bed_filename: str | Path,
    barcodes_tsv_filename: str | Path,
    blacklist_bed_filename: str | Path | None = None,
    sample_id: str | None = None,
    cb_end_to_remove: str | None = "-1",
    cb_sample_separator: str | None = "___",
):
    """
    Create fragments matrix from a fragment file and BED file with regions.
    Parameters
    ----------
    fragments_bed_filename
        Fragments BED filename.
    regions_bed_filename
        Consensus peaks / SCREEN regions BED file for which to make the fragments matrix per cell barcode.
    barcodes_tsv_filename
        TSV file with selected cell barcodes after pycisTopic QC filtering.
    blacklist_bed_filename
        BED file with blacklisted regions (Amemiya et al., 2019). Default: None.
    sample_id
        Optional sample ID to append after cell barcode after removing `cb_end_to_remove`
        and appending `cb_sample_separator`.
    cb_end_to_remove
        Remove this string from the end of the cell barcode if `sample_id` is specified.
    cb_sample_separator
        Add this string to the cell barcode if `sample_id` is specified, after removing
        `cb_end_to_remove` and before appending `sample_id`.
    Returns
    -------
    (
        counts_fragments_matrix,
        cbs,
        region_ids,
    )
    References
    ----------
    Amemiya, H. M., Kundaje, A., & Boyle, A. P. (2019). The ENCODE blacklist: identification of problematic regions of the genome. Scientific reports, 9(1), 1-5.
    """
    # Create logger
    # level = logging.INFO
    # log_format = "%(asctime)s %(name)-12s %(levelname)-8s %(message)s"
    # handlers = [logging.StreamHandler(stream=sys.stdout)]
    # logging.basicConfig(level=level, format=log_format, handlers=handlers)
    # log = logging.getLogger("cisTopic")
    # Read file with cell barcodes as a Polars Series and add sample ID to cell barcodes.
    cbs = read_barcodes_file_to_polars_series(
        barcodes_tsv_filename=barcodes_tsv_filename,
        sample_id=sample_id,
        cb_end_to_remove=cb_end_to_remove,
        cb_sample_separator=cb_sample_separator,
    )
    # log.info("Reading data for " + project)
    # Read fragments file to Polars Dataframe and add sample ID to cell barcodes.
    fragments_df_pl = read_fragments_to_polars_df(
        fragments_bed_filename=fragments_bed_filename,
        engine="pyarrow",
        sample_id=sample_id,
        cb_end_to_remove=cb_end_to_remove,
        cb_sample_separator=cb_sample_separator,
    )
    # Only keep fragments with the requested cell barcode.
    fragments_cb_filtered_df_pl = pycisTopic.fragments.filter_fragments_by_cb(
        fragments_df_pl=fragments_df_pl,
        cbs=cbs,
    ).rename({"Name": "CB", "Score": "CB_count"})
    del fragments_df_pl
    # Read regions BED file as a Polars Dataframe.
    regions_df_pl = (
        pycisTopic.fragments.read_bed_to_polars_df(
            bed_filename=regions_bed_filename,
            engine="polars",
            min_column_count=3,
        )
        .with_columns(
            (
                pl.col("Chromosome")
                + ":"
                + pl.col("Start").cast(pl.Utf8)
                + "-"
                + pl.col("End").cast(pl.Utf8)
            )
            .cast(pl.Categorical)
            .alias("RegionID")
        )
        .select(
            pl.col("Chromosome"),
            pl.col("Start"),
            pl.col("End"),
            pl.col("RegionID"),
        )
    )
    if blacklist_bed_filename:
        # Read BED file with blacklisted regions .
        blacklist_df_pl = pycisTopic.fragments.read_bed_to_polars_df(
            bed_filename=blacklist_bed_filename,
            engine="polars",
            min_column_count=3,
        ).select(
            pl.col("Chromosome"),
            pl.col("Start"),
            pl.col("End"),
        )
        # Filter out regions that overlap with blacklisted regions.
        regions_df_pl = (
            regions_df_pl.lazy()
            .join(
                # Get all regionIDs that overlap with blacklisted regions.
                pycisTopic.fragments.gr_overlap(
                    regions1_df_pl=regions_df_pl,
                    regions2_df_pl=blacklist_df_pl,
                    how="first",
                )
                .lazy()
                .select(
                    pl.col("RegionID"),
                ),
                on="RegionID",
                how="anti",
            )
            .select(
                pl.col("Chromosome"),
                pl.col("Start"),
                pl.col("End"),
                pl.col("RegionID"),
            )
            .collect()
        )
    # Get accessibility (binary and counts) for each region ID and cell barcode.
    region_cb_df_pl = (
        pycisTopic.fragments.gr_intersection(
            regions1_df_pl=regions_df_pl,
            regions2_df_pl=fragments_cb_filtered_df_pl,
            # how: Literal["all", "containment", "first", "last"] | str | None = None,
            how="all",
            regions1_info=True,
            regions2_info=True,
            regions1_coord=False,
            regions2_coord=False,
            regions1_suffix="@1",
            regions2_suffix="@2",
        )
        .rename({"CB@2": "CB"})
        .lazy()
        .group_by(["RegionID", "CB"])
        .agg(
            # Get accessibility in binary form.
            pl.lit(1).cast(pl.Int8).alias("accessible_binary"),
            # Get accessibility in count form.
            pl.len().cast(pl.UInt32).alias("accessible_count"),
        )
        .join(
            regions_df_pl.lazy()
            .select(pl.col("RegionID"))
            .with_row_index("region_idx"),
            on="RegionID",
            how="left",
        )
        .join(
            cbs.to_frame().lazy().with_row_index("CB_idx"),
            on="CB",
            how="left",
        )
        .collect()
    )
    # Construct binary accessibility matrix as a sparse matrix
    # (regions as rows and cells as columns).
    counts_fragments_matrix = sp.sparse.csr_matrix(
        (
            # All data points are 1:
            #   - same as: region_cb_df_pl.get_column("accessible_binary").to_numpy()
            #   - for count matrix: region_cb_df_pl.get_column("accessible_count").to_numpy()
            # np.ones(region_cb_df_pl.shape[0], dtype=np.int8),
            region_cb_df_pl.get_column("accessible_count").to_numpy(),
            (
                # Row indices:
                region_cb_df_pl.get_column("region_idx").to_numpy(),
                # Column indices:
                region_cb_df_pl.get_column("CB_idx").to_numpy(),
            ),
        )
    )
    return (
        counts_fragments_matrix,
        cbs.to_list(),
        regions_df_pl.get_column("RegionID").to_list(),
    )


def read_barcodes_file_to_polars_series(
    barcodes_tsv_filename: str,
    sample_id: str | None = None,
    cb_end_to_remove: str | None = "-1",
    cb_sample_separator: str | None = "___",
) -> pl.Series:
    """
    Read barcode TSV file to a Polars Series.

    Parameters
    ----------
    barcodes_tsv_filename
        TSV file with CBs.
    sample_id
        Optional sample ID to append after cell barcode after removing `cb_end_to_remove`
        and appending `cb_sample_separator`.
    cb_end_to_remove
        Remove this string from the end of the cell barcode if `sample_id` is specified.
    cb_sample_separator
        Add this string to the cell barcode if `sample_id` is specified, after removing
        `cb_end_to_remove` and before appending `sample_id`.

    Returns
    -------
    Polars Series with CBs.

    See Also
    --------
    pycisTopic.fragments.filter_fragments_by_cb

    Examples
    --------
    Read gzipped barcodes TSV file to a Polars Series.

    >>> cbs = read_barcodes_file_to_polars_series(
    ...     barcodes_tsv_filename="barcodes.tsv.gz",
    ... )

    Read uncompressed barcodes TSV file to a Polars Series.

    >>> cbs = read_barcodes_file_to_polars_series(
    ...     barcodes_tsv_filename="barcodes.tsv",
    ... )

    Read gzipped barcodes TSV file to a Polars Series and add sample ID to cell
    barcode names after removing `cb_end_to_remove` string from cell barcode and
    appending `cb_sample_separator` to the cell barcode.

    >>> cbs = read_barcodes_file_to_polars_series(
    ...     barcodes_tsv_filename="barcodes.tsv",
    ...     sample_id="sample1",
    ...     cb_end_to_remove="-1",
    ...     cb_sample_separator="___",
    ... )

    """
    cbs = (
        pl.read_csv(
            barcodes_tsv_filename,
            has_header=False,
            separator="\t",
            columns=[0],
            new_columns=["CB"],
            schema={"CB": pl.Categorical},
        )
        .filter(pl.col("CB").is_not_null())
        .unique(maintain_order=True)
    )
    # Modify cell barcode if sample ID is specified or an empty string.
    if sample_id or sample_id == "":
        separator_and_sample_id = (
            f"{cb_sample_separator + sample_id}" if cb_sample_separator else sample_id
        )
        if not cb_end_to_remove:
            # Append separator and sample ID to cell barcode.
            cbs = cbs.with_columns(
                (pl.col("CB").cast(pl.Utf8) + pl.lit(separator_and_sample_id)).cast(
                    pl.Categorical
                )
            )
        else:
            # Remove `cb_end_to_remove` from the end of the cell barcode before adding
            # separator and sample ID to cell barcode.
            cbs = cbs.with_columns(
                pl.col("CB")
                .cast(pl.Utf8)
                .str.replace(cb_end_to_remove + "$", separator_and_sample_id)
                .cast(pl.Categorical)
            )
    return cbs.to_series()


def read_fragments_to_polars_df(
    fragments_bed_filename: str,
    engine: str | Literal["polars"] | Literal["pyarrow"] = "pyarrow",
    sample_id: str | None = None,
    cb_end_to_remove: str | None = "-1",
    cb_sample_separator: str | None = "___",
) -> pl.DataFrame:
    """
    Read fragments BED file to a Polars DataFrame.

    If fragments don't have a Score column, a Score columns is created by counting
    the number of fragments with the same chromosome, start, end and CB.

    Parameters
    ----------
    fragments_bed_filename
        Fragments BED filename.
    engine
        Use Polars or pyarrow to read the fragments BED file (default: `pyarrow`).
    sample_id
        Optional sample ID to append after cell barcode after removing `cb_end_to_remove`
        and appending `cb_sample_separator`.
    cb_end_to_remove
        Remove this string from the end of the cell barcode if `sample_id` is specified.
    cb_sample_separator
        Add this string to the cell barcode if `sample_id` is specified, after removing
        `cb_end_to_remove` and before appending `sample_id`.

    Returns
    -------
    Polars DataFrame with fragments.

    See Also
    --------
    pycisTopic.fragments.read_bed_to_polars_df

    Examples
    --------
    Read gzipped fragments BED file to a Polars DataFrame.

    >>> fragments_df_pl = read_fragments_to_polars_df(
    ...     fragments_bed_filename="fragments.tsv.gz",
    ... )

    Read uncompressed fragments BED file to a Polars DataFrame.

    >>> fragments_df_pl = read_fragments_to_polars_df(
    ...     fragments_bed_filename="fragments.tsv",
    ... )

    Read gzipped fragments BED file to a Polars DataFrame and add sample ID to cell
    barcode names after removing `cb_end_to_remove` string from cell barcode and
    appending `cb_sample_separator` to the cell barcode.

    >>> fragments_df_pl = read_fragments_to_polars_df(
    ...     fragments_bed_filename="fragments.tsv.gz",
    ...     sample_id="sample1",
    ...     cb_end_to_remove="-1",
    ...     cb_sample_separator="___",
    ... )

    """
    from pycisTopic.fragments import read_bed_to_polars_df
    fragments_df_pl = read_bed_to_polars_df(
        bed_filename=fragments_bed_filename,
        engine=engine,
        min_column_count=4,
    ).lazy()
    # If no score is provided or score column is ".", generate a score column with the
    # number of fragments which have the same chromosome, start, end and CB.
    if fragments_df_pl.collect_schema().get("Score") in (None, pl.Utf8):
        fragments_df_pl = fragments_df_pl.group_by(
            ["Chromosome", "Start", "End", "Name"]
        ).agg(pl.len().cast(pl.Int32()).alias("Score"))
    else:
        fragments_df_pl = fragments_df_pl.with_columns(pl.col("Score").cast(pl.Int32()))
    # Modify cell barcode if sample ID is specified or an empty string.
    if sample_id or sample_id == "":
        separator_and_sample_id = (
            f"{cb_sample_separator + sample_id}" if cb_sample_separator else sample_id
        )
        if not cb_end_to_remove:
            # Append separator and sample ID to cell barcode.
            fragments_df_pl = fragments_df_pl.with_columns(
                (pl.col("Name").cast(pl.Utf8) + pl.lit(separator_and_sample_id)).cast(
                    pl.Categorical
                )
            )
        else:
            # Remove `cb_end_to_remove` from the end of the cell barcode before adding
            # separator and sample ID to cell barcode.
            fragments_df_pl = fragments_df_pl.with_columns(
                pl.col("Name")
                .cast(pl.Utf8)
                .str.replace(cb_end_to_remove + "$", separator_and_sample_id)
                .cast(pl.Categorical)
            )
    fragments_df_pl = fragments_df_pl.collect()
    return fragments_df_pl


def _score_regions_to_single_gene(
    acc: np.ndarray,
    exp: np.ndarray,
    gene_name: str,
    region_names: Set[str],
    regressor_type: Literal["RF", "ET", "GBM", "PR", "SR"],
    regressor_kwargs: dict,
    mask_expr_dropout: bool
    ) -> Optional[Tuple[str, pd.Series]]:
    """
    Calculates region to gene importances or region to gene correlations for a single gene

    Parameters
    ----------
    acc: np.ndarray
        Numpy array containing matrix of accessibility of regions in search space.
    exp: 
        Numpy array containing expression vector.
    gene_name: str
        Name of the gene.
    region_names: List[str]
        Names of the regions.
    regressor_type: Literal["RF", "ET", "GBM", "PR", "SR"]
        Regressor type to use, must be any of "RF", "ET", "GBM", "PR", "SR".
    regressor_kwargs: dict
        Keyword arguments to pass to the regressor function.
    mask_expr_dropout: bool
        Wether or not to mask expression dropouts.
    
    Returns
    -------
    feature_importance for regression methods and correlation_coef for correlation methods
    """
    if mask_expr_dropout:
        cell_non_zero = exp != 0
        exp = exp[cell_non_zero]
        acc = acc[cell_non_zero, :]
    # Check-up for genes with 1 region only, related to issue 2
    if acc.ndim == 1:
        acc = acc.reshape(-1, 1)
    if regressor_type in SKLEARN_REGRESSOR_FACTORY.keys():
        from arboreto import core as arboreto_core

        # fit model
        fitted_model = arboreto_core.fit_model(regressor_type=regressor_type,
                                               regressor_kwargs=regressor_kwargs,
                                               tf_matrix=acc,
                                               target_gene_expression=exp)
        # get importance scores for each feature
        feature_importance = arboreto_core.to_feature_importances(regressor_type=regressor_type,
                                                                  regressor_kwargs=regressor_kwargs,
                                                                  trained_regressor=fitted_model)
        return gene_name, pd.Series(feature_importance, index=region_names)

    elif regressor_type in SCIPY_CORRELATION_FACTORY.keys():
        # define correlation method
        correlator = SCIPY_CORRELATION_FACTORY[regressor_type]

        # do correlation and get correlation coef
        correlation_result = np.array([correlator(x, exp) for x in acc.T])
        correlation_coef = correlation_result[:, 0]

        return gene_name, pd.Series(correlation_coef, index=region_names)
    else:
        raise ValueError("Unsuported regression model")

def _get_acc_idx_per_gene(
        scplus_region_names: pd.Index,
        search_space: pd.DataFrame) -> Tuple[np.ndarray, List[List[str]]]:
    region_names = search_space["Name"].to_numpy()
    gene_names = search_space["Gene"].to_numpy()
    s = np.argsort(gene_names)
    region_names = region_names[s]
    gene_names = gene_names[s]
    region_names_to_idx = pd.DataFrame(
        index = scplus_region_names,
        data = {'idx': np.arange(len(scplus_region_names))})
    unique_gene_names, gene_idx = np.unique(gene_names, return_index = True)
    region_idx_per_gene = []
    for i in range(len(gene_idx)):
        if i < len(gene_idx) - 1:
            region_idx_per_gene.append(
                region_names_to_idx.loc[region_names[gene_idx[i]:gene_idx[i+1]], 'idx'].to_list())
        else:
            region_idx_per_gene.append(
                region_names_to_idx.loc[region_names[gene_idx[i]:], 'idx'].to_list())
    return unique_gene_names, region_idx_per_gene

def _score_regions_to_genes(
        df_exp_mtx: pd.DataFrame,
        df_acc_mtx: pd.DataFrame,
        search_space: pd.DataFrame,
        mask_expr_dropout: bool,
        regressor_type: Literal["RF", "ET", "GBM", "PR", "SR"],
        regressor_kwargs: dict,
        n_cpu: int,
        temp_dir: Union[None, pathlib.Path]) -> dict:
    """
    # TODO: Add doctstrings
    """
    if len(set(df_exp_mtx.columns)) != len(df_exp_mtx.columns):
        raise ValueError("Expression matrix contains duplicate gene names")
    if len(set(df_acc_mtx.columns)) != len(df_acc_mtx.columns):
        raise ValueError("Chromatin accessibility matrix contains duplicate gene names")
    if temp_dir is not None:
        if type(temp_dir) == str:
            temp_dir = pathlib.Path(temp_dir)
        if not temp_dir.exists():
            Warning(f"{temp_dir} does not exist, creating it.")
            os.makedirs(temp_dir)
    scplus_region_names = df_acc_mtx.columns
    scplus_gene_names = df_exp_mtx.columns
    search_space = search_space[search_space['Name'].isin(scplus_region_names)]
    search_space = search_space[search_space['Gene'].isin(scplus_gene_names)]
    # Get region indeces per gene
    gene_names, acc_idx = _get_acc_idx_per_gene(
        scplus_region_names = scplus_region_names, search_space = search_space)
    EXP = df_exp_mtx[gene_names].to_numpy()
    ACC = df_acc_mtx.to_numpy()
    regions_to_genes = dict(
        joblib.Parallel(
            n_jobs = n_cpu,
            temp_folder=temp_dir)(
                joblib.delayed(_score_regions_to_single_gene)(
                    acc = ACC[:, acc_idx[idx]],
                    exp = EXP[:, idx],
                    gene_name = gene_names[idx],
                    region_names = scplus_region_names[acc_idx[idx]],
                    regressor_type = regressor_type,
                    regressor_kwargs = regressor_kwargs, 
                    mask_expr_dropout = mask_expr_dropout
                )
                for idx in tqdm(
                    range(len(gene_names)),
                    total = len(gene_names),
                    desc=f'Running using {n_cpu} cores')
                ))
    return regions_to_genes

def calculate_regions_to_genes_relationships(
        df_exp_mtx: pd.DataFrame,
        df_acc_mtx: pd.DataFrame,
        search_space: pd.DataFrame,
        temp_dir: pathlib.Path,
        mask_expr_dropout: bool = False,
        importance_scoring_method: Literal["RF", "ET", "GBM"] = 'GBM',
        importance_scoring_kwargs: dict = GBM_KWARGS,
        correlation_scoring_method: Literal["PR", "SR"] = 'SR',
        n_cpu: int = 1,
        add_distance: bool = True):
    """
    # TODO: add docstrings
    """
    # Create logger
    level = logging.INFO
    format = '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'
    handlers = [logging.StreamHandler(stream=sys.stdout)]
    logging.basicConfig(level=level, format=format, handlers=handlers)
    log = logging.getLogger('R2G')
    # calulcate region to gene importance
    log.info(
        f'Calculating region to gene importances, using {importance_scoring_method} method')
    region_to_gene_importances = _score_regions_to_genes(
        df_exp_mtx=df_exp_mtx,
        df_acc_mtx=df_acc_mtx,
        search_space=search_space,
        mask_expr_dropout = mask_expr_dropout,
        regressor_type = importance_scoring_method,
        regressor_kwargs = importance_scoring_kwargs,
        n_cpu = n_cpu,
        temp_dir = temp_dir)

    # calculate region to gene correlation
    log.info(
        f'Calculating region to gene correlation, using {correlation_scoring_method} method')
    region_to_gene_correlation = _score_regions_to_genes(
        df_exp_mtx=df_exp_mtx,
        df_acc_mtx=df_acc_mtx,
        search_space=search_space,
        mask_expr_dropout = mask_expr_dropout,
        regressor_type = correlation_scoring_method,
        regressor_kwargs = importance_scoring_kwargs,
        n_cpu = n_cpu,
        temp_dir = temp_dir)

    # transform dictionaries to pandas dataframe
    result_df = pd.concat([pd.DataFrame(data={'target': gene,
                                                'region': region_to_gene_importances[gene].index.to_list(),
                                                'importance': region_to_gene_importances[gene].to_list(),
                                                'rho': region_to_gene_correlation[gene].loc[
                                                    region_to_gene_importances[gene].index.to_list()].to_list()})
                            for gene in region_to_gene_importances.keys()
                            ]
                            )
    result_df = result_df.reset_index()
    result_df = result_df.drop('index', axis=1)
    result_df['importance_x_rho'] = result_df['rho'] * \
        result_df['importance']
    result_df['importance_x_abs_rho'] = abs(
        result_df['rho']) * result_df['importance']
    if add_distance:
        search_space_rn = search_space.rename(
            {'Name': 'region', 'Gene': 'target'}, axis=1).copy()
        result_df = result_df.merge(search_space_rn, on=['region', 'target'])
        #result_df['Distance'] = result_df['Distance'].map(lambda x: x[0])
    log.info('Done!')
    return result_df


# changed to accept cistarget_db directly opened
def _run_cistarget_single_region_set(
        ctx_db,
        region_set,
        name,
        species,
        auc_threshold,
        nes_threshold,
        rank_threshold,
        path_to_motif_annotations,
        annotation_version,
        annotations_to_use,
        motif_similarity_fdr,
        orthologous_identity_threshold) -> pycistarget.motif_enrichment_cistarget.cisTarget:
    """Helper function to run cisTarget on a single region set."""
    cistarget_result = pycistarget.motif_enrichment_cistarget.cisTarget(
        region_set=region_set,
        name=name,
        species=species,
        auc_threshold=auc_threshold,
        nes_threshold=nes_threshold,
        rank_threshold=rank_threshold,
        path_to_motif_annotations=path_to_motif_annotations,
        annotation_version=annotation_version,
        annotation_to_use=annotations_to_use,
        motif_similarity_fdr=motif_similarity_fdr,
        orthologous_identity_threshold=orthologous_identity_threshold,
    )
    cistarget_result.run_ctx(ctx_db)
    return cistarget_result


def run_motif_enrichment_cistarget(
        region_set_dict: Dict[str, pr.PyRanges],
        ctx_db: pycistarget.motif_enrichment_cistarget.cisTargetDatabase,
        output_fname_cistarget_result: pathlib.Path,
        n_cpu: int,
        fraction_overlap_w_cistarget_database: float,
        auc_threshold: float,
        nes_threshold: float,
        rank_threshold: float,
        path_to_motif_annotations: str,
        annotation_version: str,
        motif_similarity_fdr: float,
        orthologous_identity_threshold: float,
        temp_dir: pathlib.Path,
        species: Literal[
            "homo_sapiens", "mus_musculus", "drosophila_melanogaster"],
        annotations_to_use: List[str]
):
    print("starting")
    """
    Run motif enrichment using cistarget algorithm.
    """
    cistarget_results = joblib.Parallel(
        n_jobs=n_cpu,
        temp_folder=temp_dir
    )(
        joblib.delayed(
            _run_cistarget_single_region_set
        )(
            name=key,
            region_set=region_set_dict[key],
            ctx_db=ctx_db,
            species=species,
            auc_threshold=auc_threshold,
            nes_threshold=nes_threshold,
            rank_threshold=rank_threshold,
            path_to_motif_annotations=path_to_motif_annotations,
            annotation_version=annotation_version,
            annotations_to_use=annotations_to_use,
            motif_similarity_fdr=motif_similarity_fdr,
            orthologous_identity_threshold=orthologous_identity_threshold
        )
        for key in region_set_dict
    )
    for cistarget_result in cistarget_results:
        if len(cistarget_result.motif_enrichment) > 0:
            cistarget_result.write_hdf5(
                path=output_fname_cistarget_result,
                mode="a"
            )


def run_motif_enrichment_dem(
        region_set_dict: Dict[str, pr.PyRanges],
        dem_db_fname: pathlib.Path,
        output_fname_dem_html: pathlib.Path,
        output_fname_dem_result: pathlib.Path,
        n_cpu: int,
        temp_dir: pathlib.Path,
        species: Literal[
                "homo_sapiens", "mus_musculus", "drosophila_melanogaster"],
        fraction_overlap_w_dem_database: float = 0.4,
        max_bg_regions: Optional[int] = None,
        path_to_genome_annotation: Optional[str] = None,
        balance_number_of_promoters: bool = True,
        promoter_space: int = 1_000,
        adjpval_thr: float = 0.05,
        log2fc_thr: float = 1.0,
        mean_fg_thr: float = 0.0,
        motif_hit_thr: Optional[float] = None,
        path_to_motif_annotations: Optional[str] = None,
        annotation_version: str = "v10nr_clust",
        annotations_to_use: tuple = ("Direct_annot", "Orthology_annot"),
        motif_similarity_fdr: float = 0.001,
        orthologous_identity_threshold: float = 0.0,
        seed: int = 555,
        write_html: bool = True):
    """
    Run motif enrichment using DEM algorithm.

    region_set_folder --> replaced to region_set_dictionary
    """
    from pycistarget.motif_enrichment_dem import (
        DEM,
    )

    # Read genome annotation, if needed
    if path_to_genome_annotation is not None:
        genome_annotation = pd.read_table(path_to_genome_annotation)
    else:
        genome_annotation = None

    dem_results: List[DEM] = joblib.Parallel(
        n_jobs=n_cpu,
        temp_folder=temp_dir
    )(
        joblib.delayed(
            _run_dem_single_region_set
        )(
            foreground_region_sets=foreground_region_sets,
            background_region_sets=background_region_sets,
            name=name,
            dem_db_fname=dem_db_fname,
            max_bg_regions=max_bg_regions,
            genome_annotation=genome_annotation,
            balance_number_of_promoters=balance_number_of_promoters,
            promoter_space=promoter_space,
            seed=seed,
            fraction_overlap_w_dem_database=fraction_overlap_w_dem_database,
            species=species,
            adjpval_thr=adjpval_thr,
            log2fc_thr=log2fc_thr,
            mean_fg_thr=mean_fg_thr,
            motif_hit_thr=motif_hit_thr,
            path_to_motif_annotations=path_to_motif_annotations,
            annotation_version=annotation_version,
            annotations_to_use=annotations_to_use,
            motif_similarity_fdr=motif_similarity_fdr,
            orthologous_identity_threshold=orthologous_identity_threshold
        )
        for name, foreground_region_sets, background_region_sets in _get_foreground_background(region_set_dict)
    )
    if write_html:
        all_motif_enrichment_df = pd.concat(
            ctx_result.motif_enrichment for ctx_result in dem_results
        )
        all_motif_enrichment_df.to_html(
            buf = output_fname_dem_html,
            escape = False,
            col_space = 80
        )
    for dem_result in dem_results:
        if len(dem_result.motif_enrichment) > 0:
            print("Saving cistrome.")
            dem_result.write_hdf5(
                path = output_fname_dem_result,
                mode = "a"
            )
        else:
            warnings.warn("Warning.....................No cistrome selected, you should probably lower dem_adj_pval_thr if you wanna keep cistromes")


def _run_dem_single_region_set(
        foreground_region_sets,
        background_region_sets,
        dem_db_fname,
        max_bg_regions,
        genome_annotation,
        balance_number_of_promoters,
        promoter_space,
        seed,
        fraction_overlap_w_dem_database,
        name,
        species,
        adjpval_thr,
        log2fc_thr,
        mean_fg_thr,
        motif_hit_thr,
        path_to_motif_annotations,
        annotation_version,
        annotations_to_use,
        motif_similarity_fdr,
        orthologous_identity_threshold):  # -> DEM:
    """Helper function to run DEM on a single region set."""
    from pycistarget.motif_enrichment_dem import (
        DEMDatabase,
        get_foreground_and_background_regions,
    )
    # Get foreground and background regions for DEM analysis
    foreground_regions, background_regions = get_foreground_and_background_regions(
        foreground_region_sets = foreground_region_sets,
        background_region_sets = background_region_sets,
        max_bg_regions = max_bg_regions,
        genome_annotation = genome_annotation,
        balance_number_of_promoters = balance_number_of_promoters,
        promoter_space = promoter_space,
        seed = seed)

    # Load DEM database
    dem_db = DEMDatabase(
        dem_db_fname,
        fraction_overlap=fraction_overlap_w_dem_database)
    # Setup DEM analysis
    dem_result = DEM(
        foreground_regions = foreground_regions,
        background_regions = background_regions,
        name = name,
        species = species,
        adjpval_thr = adjpval_thr,
        log2fc_thr = log2fc_thr,
        mean_fg_thr = mean_fg_thr,
        motif_hit_thr = motif_hit_thr,
        path_to_motif_annotations = path_to_motif_annotations,
        annotation_version = annotation_version,
        annotation_to_use = annotations_to_use,
        motif_similarity_fdr = motif_similarity_fdr,
        orthologous_identity_threshold = orthologous_identity_threshold)
    # Run DEM analysis
    dem_result.run(dem_db)
    return dem_result


def prepare_motif_enrichment_results(
        paths_to_motif_enrichment_results: List[str],
        mdata: mudata.MuData,
        out_file_direct_annotation: pathlib.Path,
        out_file_extended_annotation: pathlib.Path,
        out_file_tf_names: pathlib.Path,
        direct_annotation: List[str],
        extended_annotation: List[str]) -> None:
    """
    Prepare motif enrichment results for SCENIC+ analysis.

    Parameters
    ----------
    paths_to_motif_enrichment_results : List[str]
        List of paths to motif enrichment results.
    multiome_mudata_fname : pathlib.Path
        Path to multiome MuData file.
    out_file_direct_annotation : pathlib.Path
        Path to store TF cistromes with direct annotation.
    out_file_extended_annotation : pathlib.Path
        Path to store TF cistromes with extended annotation.
    out_file_tf_names : pathlib.Path
        Path to store TF names.
    direct_annotation : List[str]
        List of annotations to use for direct annotation.
    extended_annotation : List[str]
        List of annotations to use for extended annotation.

    """
    from scenicplus.data_wrangling.cistarget_wrangling import get_and_merge_cistromes
    log.info("Reading multiome MuData.")
    log.info("Getting cistromes.")
    region_names = mdata["atac"].var_names.str.replace('-', ':', n=1)  # replace first '-' by ':'

    adata_direct_cistromes, adata_extended_cistromes = get_and_merge_cistromes(
        paths_to_motif_enrichment_results=paths_to_motif_enrichment_results,
        scplus_regions=set(region_names),
        direct_annotation=direct_annotation,
        extended_annotation=extended_annotation)
    # Get transcription factor names from cistromes
    # Later, to calculate TF-to-gene relationships these TFs will be used.
    TFs = list(
        {
            *adata_direct_cistromes.var_names,
            *adata_extended_cistromes.var_names
        } & set(mdata["rna"].var_names))
    log.info(f"Found {len(TFs)} TFs.")
    log.info(f"Saving TF names to: {out_file_tf_names.__str__()}")
    with open(out_file_tf_names, "w") as f:
        for TF in TFs:
            _ = f.write(TF)
            _ = f.write("\n")
    if len(direct_annotation) > 0:
        log.info(
            f"Writing direct cistromes to: {out_file_direct_annotation.__str__()}")
        adata_direct_cistromes.write_h5ad(out_file_direct_annotation.__str__())
    if len(extended_annotation) > 0:
        log.info(
            f"Writing extended cistromes to: {out_file_extended_annotation.__str__()}")
        adata_extended_cistromes.write_h5ad(out_file_extended_annotation.__str__())



def infer_TF_to_gene(
        multiome_mudata_fname: pathlib.Path,
        tf_names: pathlib.Path,
        temp_dir: pathlib.Path,
        method: Literal["GBM", "RF"],
        n_cpu: int,
        seed: int):
    """
    Replace tf_names_fname by tf_list 
    Replace scRNA by rna
    """
    from scenicplus.TF_to_gene import calculate_TFs_to_genes_relationships
    log.info("Reading multiome MuData.")
    mdata = mudata.read(multiome_mudata_fname.__str__())
    log.info(f"Using {len(tf_names)} TFs.")
    adj = calculate_TFs_to_genes_relationships(
        df_exp_mtx=mdata["rna"].to_df(),
        tf_names = tf_names,
        temp_dir = temp_dir,
        method = method,
        n_cpu = n_cpu,
        seed = seed)

    return adj


def infer_grn(
        TF_to_gene_adj: pathlib.Path,
        region_to_gene_adj: pathlib.Path,
        cistromes: pathlib.Path,
        ranking_db_fname: str,
        is_extended: bool,
        temp_dir: pathlib.Path,
        order_regions_to_genes_by: str,
        order_TFs_to_genes_by: str,
        gsea_n_perm: int,
        quantiles: List[float],
        top_n_regionTogenes_per_gene: List[float],
        top_n_regionTogenes_per_region: List[float],
        binarize_using_basc: bool,
        min_regions_per_gene: int,
        rho_dichotomize_tf2g: bool,
        rho_dichotomize_r2g: bool,
        rho_dichotomize_eregulon: bool,
        keep_only_activating: bool,
        rho_threshold: float,
        min_target_genes: int,
        n_cpu: int):
    """
    Pass TF_to_gene_adj, region_to_gene_adj, cistromes directly
    """
    #from scenicplus.triplet_score import calculate_triplet_score
    log.info("Loading TF to gene adjacencies.")
    tf_to_gene = TF_to_gene_adj

    log.info("Loading region to gene adjacencies.")
    region_to_gene = region_to_gene_adj

    log.info("Loading cistromes.")
    cistromes = cistromes

    eRegulons = build_grn(
        tf_to_gene=tf_to_gene,
        region_to_gene=region_to_gene,
        cistromes=cistromes,
        is_extended=is_extended,
        temp_dir=temp_dir.__str__(),
        order_regions_to_genes_by=order_regions_to_genes_by,
        order_TFs_to_genes_by=order_TFs_to_genes_by,
        gsea_n_perm=gsea_n_perm,
        quantiles=quantiles,
        top_n_regionTogenes_per_gene=top_n_regionTogenes_per_gene,
        top_n_regionTogenes_per_region=top_n_regionTogenes_per_region,
        binarize_using_basc=binarize_using_basc,
        min_regions_per_gene=min_regions_per_gene,
        rho_dichotomize_tf2g=rho_dichotomize_tf2g,
        rho_dichotomize_r2g=rho_dichotomize_r2g,
        rho_dichotomize_eregulon=rho_dichotomize_eregulon,
        keep_only_activating=keep_only_activating,
        rho_threshold=rho_threshold,
        NES_thr=0,
        adj_pval_thr=1,
        min_target_genes=min_target_genes,
        n_cpu=n_cpu,
        merge_eRegulons=True,
        disable_tqdm=False)

    log.info("Formatting eGRN as table.")
    eRegulon_metadata = _format_egrns(
        eRegulons=eRegulons,
        tf_to_gene=tf_to_gene)

    print(eRegulon_metadata)
    log.info("Calculating triplet ranking.")
    #print(TF_to_region_score)
   # eRegulon_metadata = calculate_triplet_score(
    #    TF_to_region_score=TF_to_region_score,
     #   eRegulon_metadata=eRegulon_metadata)

    #log.info(f"Saving network to {eRegulon_out_fname.__str__()}")
    return eRegulon_metadata


def build_grn(
        tf_to_gene: pd.DataFrame,
        region_to_gene: pd.DataFrame,
        cistromes: anndata.AnnData,
        is_extended: bool,
        temp_dir: str,
        order_regions_to_genes_by='importance',
        order_TFs_to_genes_by='importance',
        gsea_n_perm=1000,
        quantiles=(0.85, 0.90),
        top_n_regionTogenes_per_gene=(5, 10, 15),
        top_n_regionTogenes_per_region=(),
        binarize_using_basc=True,
        min_regions_per_gene=0,
        rho_dichotomize_tf2g=True,
        rho_dichotomize_r2g=True,
        rho_dichotomize_eregulon=True,
        keep_only_activating=False,
        rho_threshold=RHO_THRESHOLD,
        NES_thr=0,
        adj_pval_thr=1,
        min_target_genes=5,
        n_cpu=1,
        merge_eRegulons=True,
        disable_tqdm=False,
        **kwargs) -> List[eRegulon]:
    log.info('Thresholding region to gene relationships')
    # some tfs are missing from tf_to_gene because they are not 
    # preset in the gene expression matrix, so subset!
    cistromes = cistromes[
        :, cistromes.var_names[cistromes.var_names.isin(tf_to_gene['TF'])]]
    relevant_tfs, e_modules = create_emodules(
        region_to_gene=region_to_gene,
        cistromes=cistromes,
        is_extended=is_extended,
        order_regions_to_genes_by=order_regions_to_genes_by,
        quantiles=quantiles,
        top_n_regionTogenes_per_gene=top_n_regionTogenes_per_gene,
        top_n_regionTogenes_per_region=top_n_regionTogenes_per_region,
        binarize_using_basc=binarize_using_basc,
        min_regions_per_gene=min_regions_per_gene,
        rho_dichotomize=rho_dichotomize_r2g,
        keep_only_activating=keep_only_activating,
        rho_threshold=rho_threshold,
        disable_tqdm=disable_tqdm,
        n_cpu=n_cpu,
        temp_dir=temp_dir)
    log.info('Subsetting TF2G adjacencies for TF with motif.')
    TF2G_adj_relevant = tf_to_gene.loc[tf_to_gene['TF'].isin(relevant_tfs)]
    TF2G_adj_relevant.index = TF2G_adj_relevant["TF"]
    log.info('Running GSEA...')
    if rho_dichotomize_tf2g:
        log.info("Generating rankings...")
        TF2G_adj_relevant_pos = TF2G_adj_relevant.loc[TF2G_adj_relevant["rho"] > rho_threshold]
        TF2G_adj_relevant_neg = TF2G_adj_relevant.loc[TF2G_adj_relevant["rho"] < -rho_threshold]
        pos_TFs, c = np.unique(TF2G_adj_relevant_pos["TF"], return_counts=True)
        pos_TFs = pos_TFs[c >= min_target_genes]
        neg_TFs, c = np.unique(TF2G_adj_relevant_neg["TF"], return_counts=True)
        neg_TFs = neg_TFs[c >= min_target_genes]
        print(len(pos_TFs), len(neg_TFs))
        # The expression below will fail if there is only a single target gene (after thresholding on rho)
        # TF2G_adj_relevant_pos/neg.loc[TF] will return a pd.Series instead of dataframe
        # This should never be the case though (if min_target_genes > 1)
        # But better fix this at some point!
        TF_to_ranking_pos = {
            TF: TF2G_adj_relevant_pos.loc[TF].set_index('target')[order_TFs_to_genes_by].sort_values(ascending = False)
            for TF in tqdm(pos_TFs, total = len(pos_TFs))}
        TF_to_ranking_neg = {
            TF: TF2G_adj_relevant_neg.loc[TF].set_index('target')[order_TFs_to_genes_by].sort_values(ascending = False)
            for TF in tqdm(neg_TFs, total = len(neg_TFs))}
        pos_tf_gene_modules = joblib.Parallel(
            n_jobs=n_cpu,
            temp_folder=temp_dir)(
            joblib.delayed(_run_gsea_for_e_module)(
                e_module=e_module,
                rnk=TF_to_ranking_pos[e_module.transcription_factor],
                gsea_n_perm=gsea_n_perm,
                context=frozenset(['positive tf2g']))
            for e_module in tqdm(
                e_modules, 
                total = len(e_modules),
                desc="Running for Positive TF to gene")
            if e_module.transcription_factor in pos_TFs)
        neg_tf_gene_modules = joblib.Parallel(
            n_jobs=n_cpu,
            temp_folder=temp_dir)(
            joblib.delayed(_run_gsea_for_e_module)(
                e_module=e_module,
                rnk=TF_to_ranking_neg[e_module.transcription_factor],
                gsea_n_perm=gsea_n_perm,
                context=frozenset(['negative tf2g']))
            for e_module in tqdm(
                e_modules, 
                total = len(e_modules),
                desc="Running for Negative TF to gene")
            if e_module.transcription_factor in neg_TFs)
        new_e_modules = [*pos_tf_gene_modules, *neg_tf_gene_modules]
    else:
        log.info("Generating rankings...")
        TFs, c = np.unique(TF2G_adj_relevant["TF"], return_counts=True)
        TFs = TFs[c >= min_target_genes]
        # The expression below will fail if there is only a single target gene (after thresholding on rho)
        # TF2G_adj_relevant.loc[TF] will return a pd.Series instead of dataframe
        # This should never be the case though (if min_target_genes > 1)
        # But better fix this at some point!
        TF_to_ranking = {
            TF: TF2G_adj_relevant.loc[TF].set_index('target')[order_TFs_to_genes_by].sort_values(ascending = False)
            for TF in tqdm(TFs, total = len(TFs))}
        new_e_modules = joblib.Parallel(
            n_jobs=n_cpu,
            temp_folder=temp_dir)(
            joblib.delayed(_run_gsea_for_e_module)(
                e_module=e_module,
                rnk=TF_to_ranking[e_module.transcription_factor],
                gsea_n_perm=gsea_n_perm,
                context=frozenset(['negative tf2g']))
            for e_module in tqdm(
                e_modules, 
                total = len(e_modules),
                desc="Running for Negative TF to gene")
            if e_module.transcription_factor in TFs)
    # filter out nans
    new_e_modules = [m for m in new_e_modules if not np.isnan(
        m.gsea_enrichment_score) and not np.isnan(m.gsea_pval)]

    log.info(
        f'Subsetting on adjusted pvalue: {adj_pval_thr}, minimal NES: {NES_thr} and minimal leading edge genes {min_target_genes}')
    # subset on adj_p_val
    adj_pval = p_adjust_bh([m.gsea_pval for m in new_e_modules])
    if any([np.isnan(p) for p in adj_pval]):
        Warning(
            'Something went wrong with calculating adjusted p values, early returning!')
        return new_e_modules

    for module, adj_pval in zip(new_e_modules, adj_pval):
        module.gsea_adj_pval = adj_pval

    e_modules_to_return: List[eRegulon] = []
    for module in new_e_modules:
        if module.gsea_adj_pval < adj_pval_thr and module.gsea_enrichment_score > NES_thr:
            module_in_LE = module.subset_leading_edge(inplace=False)
            if module_in_LE.n_target_genes >= min_target_genes:
                e_modules_to_return.append(module_in_LE)
    if merge_eRegulons:
        log.info('Merging eRegulons')
        e_modules_to_return = merge_emodules(
            e_modules=e_modules_to_return, inplace=False, rho_dichotomize=rho_dichotomize_eregulon)
    e_modules_to_return = [
        x for x in e_modules_to_return if not isinstance(x, list)]
    return e_modules_to_return


def calculate_triplet_score(
        TF_to_region_score: pd.DataFrame,
        eRegulon_metadata: pd.DataFrame,
        ) -> pd.DataFrame:
    """
        Calculate the triplet score for each eRegulon.
        Adapted from SCENIC+ to take tfb scores instead of rank comigng from cistarget_db.
        TF_to_region_score: pd.DataFrame, should contain the following columns: 
                - tf
                - region
                - score
    """

    eRegulon_metadata = eRegulon_metadata.copy()

    TF_region_iter = eRegulon_metadata[["TF", "Region"]].to_numpy()
    eRegulon_mask = [(TF_to_region_score["tf"] == TF) &  (TF_to_region_score["region"] == region) for TF, region in TF_region_iter]
    eRegulon_metadata = eRegulon_metadata.iloc[eRegulon_mask, :]
    print(eRegulon_metadata)

    TF_region_iter = eRegulon_metadata[["TF", "Region"]].to_numpy()
    TF_to_region_score = np.array([
        TF_to_region_score[
            (TF_to_region_score["tf"] == TF) &
            (TF_to_region_score["region"] == region)]["score"].values[0]
        for TF, region in TF_region_iter])

    TF_to_gene_score = eRegulon_metadata["importance_TF2G"].to_numpy()
    region_to_gene_score = eRegulon_metadata["importance_R2G"].to_numpy()
    # rank the scores
    TF_to_region_rank = _rank_scores_and_assign_random_ranking_in_range_for_ties(
        TF_to_region_score) # positive since higher score is better here
    TF_to_gene_rank = _rank_scores_and_assign_random_ranking_in_range_for_ties(
        TF_to_gene_score)
    region_to_gene_rank = _rank_scores_and_assign_random_ranking_in_range_for_ties(
        region_to_gene_score)
    # create rank ratios
    TF_to_gene_rank_ratio = (TF_to_gene_rank.astype(np.float64) + 1) / TF_to_gene_rank.shape[0]
    region_to_gene_rank_ratio = (region_to_gene_rank.astype(np.float64) + 1) / region_to_gene_rank.shape[0]
    TF_to_region_rank_ratio = (TF_to_region_rank.astype(np.float64) + 1) / TF_to_region_rank.shape[0]
    # create aggregated rank
    rank_ratios = np.array([
        TF_to_gene_rank_ratio, region_to_gene_rank_ratio, TF_to_region_rank_ratio])
    aggregated_rank = np.zeros((rank_ratios.shape[1],), dtype = np.float64)
    for i in range(rank_ratios.shape[1]):
        aggregated_rank[i] = _calculate_cross_species_rank_ratio_with_order_statistics(rank_ratios[:, i])
    eRegulon_metadata["triplet_rank"] = aggregated_rank.argsort().argsort()
    return eRegulon_metadata




# Init args
parser = argparse.ArgumentParser()
    # Data and folders
parser.add_argument('-f', '--frags',  nargs='+', help='Path to fragments file')
parser.add_argument('-i', '--mudata', type=str, help='Path to metadata file')
parser.add_argument('-o', '--out', type=str, help='Path to output directory')
parser.add_argument('-t', '--temp_dir', type=str, help='Path to temp directory')
parser.add_argument('--ray_tmp_dir', type=str, help='Path to ray tmp directory')
parser.add_argument('--tmp_scenicplus', type=str, help='Path to tmp scenicplus directory')
    # tfb intermediate files
parser.add_argument('--annotation_direct_path', type=str, help='Path to direct annotation')
parser.add_argument('--annotation_extended_path', type=str, help='Path to extended annotation')
parser.add_argument('--tf_names_path', type=str, help='Path to TF names')
parser.add_argument('--cistarget_results_path', required=True)
parser.add_argument('--dem_results_path', required=True)
    # General parameters
parser.add_argument('-g', '--organism', type=str, help='Organism')
parser.add_argument('-c', '--njobs', type=int, help='Number of cores')
    # Additional files
parser.add_argument('-m', '--chrom_sizes_m', type=str, help='Path to mouse chrom sizes')
parser.add_argument('-j', '--chrom_sizes_h', type=str, help='Path to human chrom sizes')
parser.add_argument('-a', '--annot_mouse', type=str, help='Path to mouse annotations')
parser.add_argument('-b', '--annot_human', type=str, help='Path to human annotations')
parser.add_argument('-r', '--cistarget_rankings_human', type=str, help='Path to human cistarget rankings')
parser.add_argument('-s', '--cistarget_scores_human', type=str, help='Path to human cistarget scores')
parser.add_argument('--path_to_motif_annotations_human', type=str, help='Path to human motif annotations')
parser.add_argument('--path_to_motif_annotations_mouse', type=str, help='Path to mouse motif annotations')
    # ScenicPlus parameters
parser.add_argument('--search_space_upstream', type=str, help='Search space upstream')
parser.add_argument('--search_space_downstream', type=str, help='Search space downstream')
parser.add_argument('--search_space_extend_tss', type=str, help='Search space extend tss')
parser.add_argument('--remove_promoters', type=bool, help='Remove promoters')
parser.add_argument('--use_gene_boundaries', type=bool, help='Use gene boundaries')
parser.add_argument('--region_to_gene_importance_method', type=str, help='Region to gene importance method')
parser.add_argument('--region_to_gene_correlation_method', type=str, help='Region to gene correlation method')
parser.add_argument('--method_mdl', type=str, help='Method mdl')
parser.add_argument('--order_regions_to_genes_by', type=str, required=False, default="importance")
parser.add_argument('--order_TFs_to_genes_by', type=str, required=False, default="importance")
parser.add_argument('--gsea_n_perm', type=int, required=False, default=1000)
parser.add_argument('--quantile_thresholds_region_to_gene', type=float, required=False, nargs="*", default=[0.85, 0.90, 0.95])
parser.add_argument('--top_n_regionTogenes_per_gene', type=int, required=False, nargs="*", default=[5, 10, 15])
parser.add_argument('--top_n_regionTogenes_per_region', type=int, required=False, nargs="*", default=[])
parser.add_argument('--min_regions_per_gene', required=True, type=int)
parser.add_argument('--rho_threshold', type=float, required=False, default=0.05,)
parser.add_argument('--min_target_genes', type=int, required=False, default=10,)
args = vars(parser.parse_args())

# Set up variables
frags = args['frags']
mudata_file = args['mudata']
output_fname = args['out']
tmp_scenicplus = args['tmp_scenicplus']
ray_tmp_dir = os.path.join(args['ray_tmp_dir'])
njobs = args['njobs']
organism = args['organism']
dem_results_path = args['dem_results_path']
cistarget_results_path = args['cistarget_results_path']

if organism == 'hg38':
    organism = "hsapiens"
    species = "homo_sapiens",
    chromsizes_fname = args['chrom_sizes_h']
    annot_fname = args['annot_human']
    cistarget_ranking_db_fname = pathlib.Path(args["cistarget_rankings_human"])
    cistarget_score_db_fname = pathlib.Path(args["cistarget_scores_human"])
    path_to_motif_annotations = args["path_to_motif_annotations_human"]
elif organism == 'mm10':
    organism = "mmusculus"
    chromsizes_fname = args['chrom_sizes_m']
    annot_fname = args['annot_mouse']
    cistarget_ranking_db_fname = pathlib.Path(args["cistarget_rankings_mouse"])
    cistarget_score_db_fname = pathlib.Path(args["cistarget_scores_mouse"])
    path_to_motif_annotations = args["path_to_motif_annotations_mouse"]
    species = "mus_musculus"
annotation_version = "v10nr_clust"

annot = pd.read_csv(annot_fname, sep="\t")
annot.Start = annot.Start.astype(np.int32)
annot.Gene = annot.Gene.astype(str)
annot = pl.DataFrame(annot)
print(annot.head())


# cisTopic paths
cis_topic_tmp_dir = os.path.join(tmp_scenicplus, 'cisTopic')
quality_control_dir = os.path.join(cis_topic_tmp_dir, 'quality_control')
bed_folder = os.path.join(cis_topic_tmp_dir, "pseudobulk_bed_files")
bigwig_folder = os.path.join(cis_topic_tmp_dir, "pseudobulk_bw_files")
bed_pickle = os.path.join(cis_topic_tmp_dir, "pseudobulk_bw_files", 'bed_paths.pkl')
bw_pickle = os.path.join(cis_topic_tmp_dir, "pseudobulk_bw_files", 'bw_paths.pkl')
# MACS paths
macs_folder = os.path.join(ray_tmp_dir, "MACS")
narrow_peaks_pickle = os.path.join(macs_folder, 'narrow_peaks_dict.pkl')
# Consensus peaks paths
consensus_peaks_bed = os.path.join(cis_topic_tmp_dir, 'consensus_region.bed')
quality_control_dir = os.path.join(cis_topic_tmp_dir, 'quality_control')
# cisTopic object path
cistopic_obj_path = os.path.join(cis_topic_tmp_dir, 'cistopic_obj.pkl')

os.makedirs(tmp_scenicplus, exist_ok=True)
os.makedirs(cis_topic_tmp_dir, exist_ok=True)
os.makedirs(quality_control_dir, exist_ok=True)
os.makedirs(macs_folder, exist_ok=True)
os.makedirs(bed_folder, exist_ok=True)
os.makedirs(bigwig_folder, exist_ok=True)

# Celltype annotations from MuData object
mdata = mu.read_h5mu(mudata_file)
cell_data = mdata.obs
cell_data['celltype'] = cell_data['celltype'].astype(str) 
cell_data['celltype'] = cell_data['celltype'].str.replace(' ', '_')

# Create fragments dictionary
fragments_dict = {}
batch_ids = cell_data['batch'].unique()
cell_data['barcode'] = cell_data.index
fragments_files = {frag.split(".frags")[0].split("/")[-1]:frag for frag in frags}
print("fragments_files:{}".format(fragments_files))

for batch_id in batch_ids:
    # Add batch_id to fragments_dict
    fragments_dict[batch_id] = fragments_files[batch_id]
    # create tabix file
    print("Create index files ('{fragments_file}.tbi')")
    try:
        pysam.tabix_index(fragments_files[batch_id], preset="bed")
    except OSError:
        print("It seems there is already a file called '{}.tbi'. Skipping index file creation.".format(
            fragments_files[batch_id]))
    # Remove batch_id from cell barcodes
    cell_data.loc[cell_data['batch'] == batch_id, 'barcode'] = \
        cell_data.loc[cell_data['batch'] == batch_id, 'barcode'].str.replace(
            batch_id + '_', '')

cell_data['barcode'] = cell_data['barcode'].str.split('-1').str[0] + '-1'
cell_data.index = cell_data['barcode']


# Get chromosome sizes (for hg38 here)
chromsizes = pd.read_csv(chromsizes_fname, sep='\t', header=None)
chromsizes.columns = ['Chromosome', 'End']
chromsizes['Start'] = [0]*chromsizes.shape[0]
chromsizes = chromsizes.loc[:, ['Chromosome', 'Start', 'End']]
chromsizes['Chromosome'] = [chromsizes['Chromosome'][x].replace('v', '.')
                            for x in range(len(chromsizes['Chromosome']))]
chromsizes['Chromosome'] = [chromsizes['Chromosome'][x].split('_')[1]
                            if len(chromsizes['Chromosome'][x].split('_')) > 1
                            else chromsizes['Chromosome'][x]
                            for x in range(len(chromsizes['Chromosome']))]
chromsizes = pr.PyRanges(chromsizes)

bw_paths, bed_paths = pycisTopic.pseudobulk_peak_calling.export_pseudobulk(
    input_data=cell_data,
    variable='celltype',               # variable by which to generate pseubulk profiles, in this case we want pseudobulks per celltype
    sample_id_col='batch',
    chromsizes=chromsizes,
    bed_path=bed_folder,               # specify where pseudobulk_bed_files should be stored
    bigwig_path=bigwig_folder,         # specify where pseudobulk_bw_files should be stored
    path_to_fragments=fragments_dict,  # location of fragment fiels
    n_cpu=1,                           # specify the number of cores to use, we use ray for multi processing
    normalize_bigwig=True,
    temp_dir=ray_tmp_dir,              # specify the location of the temp directory
    split_pattern='-')

# Run peak calling
macs_path = 'macs2'
narrow_peaks_dict = pycisTopic.pseudobulk_peak_calling.peak_calling(
    macs_path,
    bed_paths,
    macs_folder,
    genome_size='hs',
    n_cpu=njobs,
    input_format='BEDPE',
    shift=73,
    ext_size=146,
    keep_dup='all',
    q_value=0.05,
    _temp_dir=ray_tmp_dir
    )


peak_half_width = 250
# Get consensus peaks
consensus_peaks = get_consensus_peaks(
    narrow_peaks_dict,
    peak_half_width,
    chromsizes=chromsizes,
    path_to_blacklist=None  # os.path.join(folder_path_to_blacklist, 'peak-blacklist.v2.bed')
)

consensus_peaks.to_bed(
    path=consensus_peaks_bed,
    keep=True,
    compression='infer',
    chain=False)

# transform to pl.DataFrame
consensus_peaks = consensus_peaks.as_df()
consensus_peaks = pl.DataFrame(consensus_peaks)

##################
# Get annotation #
##################
annot = pd.read_csv(annot_fname, sep='\t')
annot.Start = annot.Start.astype(np.int32)
annot.Gene = annot.Gene.astype(str)
annot = pl.DataFrame(annot)
print(annot)
print(chromsizes)

# Run QC
for key in fragments_dict.keys():
    qc(
        fragments_tsv_filename=fragments_dict[key],
        regions_bed_filename=consensus_peaks_bed,
        tss_annotation_bed_df_pl=annot,
        output_prefix=os.path.join(quality_control_dir, key),
        tss_flank_window=2000,
        tss_smoothing_rolling_window=10,
        tss_minimum_signal_window=100,
        tss_window=50,
        tss_min_norm=0.2,
        use_genomic_ranges=True,
        min_fragments_per_cb=10,
        collapse_duplicates=True,
        no_threads=njobs,
        engine='pyarrow'
    )

# Get barcodes passing QC
sample_id_to_barcodes_passing_filters = {}
sample_id_to_thresholds = {}
for sample_id in fragments_dict:
    (
        sample_id_to_barcodes_passing_filters[sample_id],
        sample_id_to_thresholds[sample_id]
    ) = get_barcodes_passing_qc_for_sample(
            sample_id=sample_id,
            pycistopic_qc_output_dir=quality_control_dir,
            unique_fragments_threshold=None,  # use automatic thresholding
            tss_enrichment_threshold=None,  # use automatic thresholding
            frip_threshold=0,
            use_automatic_thresholds=True,
    )
    pd.Series(sample_id_to_barcodes_passing_filters[sample_id]).to_csv(
        os.path.join(quality_control_dir, f'{sample_id}.barcodes_passing_qc.txt'),
        header=None, index=None, sep='\t'
    )

fragments_df_pl = pycisTopic.fragments.read_bed_to_polars_df(
        bed_filename=fragments_dict[sample_id],
        engine="polars",
        min_column_count=4,
    ).lazy()


# Create cistopic object
cistopic_obj_list = []
for sample_id in fragments_dict:
    # compute matrix from fragments, peaks, and barcodes passing QC.
    counts_fragments_matrix, cbs, region_ids = create_fragment_matrix_from_fragments(
        fragments_bed_filename=fragments_dict[sample_id],
        regions_bed_filename=consensus_peaks_bed,
        barcodes_tsv_filename=os.path.join(quality_control_dir, f'{sample_id}.barcodes_passing_qc.txt'),
    )
    cistopic_obj = pycisTopic.cistopic_class.create_cistopic_object(
        fragment_matrix=counts_fragments_matrix,
        cell_names=cbs,
        region_names=region_ids,
        path_to_fragments=fragments_dict[sample_id],
        project=sample_id
    )
    cistopic_obj_list.append(cistopic_obj)

cistopic_obj = pycisTopic.cistopic_class.merge(cistopic_obj_list)
cistopic_obj.add_cell_data(cell_data, split_pattern='-')

print("Generate models.")

models = pycisTopic.lda_models.run_cgs_models(
    cistopic_obj,
    n_topics=[40,50], #[2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50],
    n_cpu=njobs,
    n_iter=15, #150,
    random_state=555,
    alpha=50,
    alpha_by_topic=True,
    eta=0.1,
    eta_by_topic=False,
    save_path=cis_topic_tmp_dir,
    _temp_dir=ray_tmp_dir
)

model = pycisTopic.lda_models.evaluate_models(
    models,
    select_model=None,
    return_model=True
)
cistopic_obj.add_LDA_model(model)

# Impute accessibility
imputed_acc_obj = impute_accessibility(
    cistopic_obj,
    selected_cells=None,
    selected_regions=None,
    scale_factor=10**6
)
normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)

# Get top markers
variable_regions = find_highly_variable_features(
    normalized_imputed_acc_obj,
    min_disp=0.05,
    min_mean=0.0125,
    max_mean=3,
    max_disp=np.inf,
    n_bins=20,
    n_top_features=None,
    plot=True
)
markers_dict = find_diff_features(
    cistopic_obj,
    imputed_acc_obj,
    variable='celltype',
    var_features=variable_regions,
    contrasts=None,
    adjpval_thr=0.05,
    log2fc_thr=np.log2(1.5),
    n_cpu=5,
    _temp_dir=ray_tmp_dir,
    split_pattern='-'
)
region_bin_topics_otsu = binarize_topics(
    cistopic_obj, method='otsu',
    plot=True, num_columns=5
)
region_bin_topics_top_3k = binarize_topics(
    cistopic_obj, method='ntop', ntop=3_000,
    plot=True, num_columns=5
)


region_sets = {}
for topic in markers_dict:
    region_sets[topic + "markers_dict"] = \
        pr.PyRanges(
            region_names_to_coordinates(
                markers_dict[topic].index
            ).sort_values(
                ["Chromosome", "Start", "End"]
            )
        )
for topic in region_bin_topics_otsu:
    region_sets[topic + "_otsu"] = \
        pr.PyRanges(
            region_names_to_coordinates(
                region_bin_topics_otsu[topic].index
            ).sort_values(
                ["Chromosome", "Start", "End"]
            )
        )
for topic in region_bin_topics_top_3k:
    region_sets[topic + "_top3k"] = \
        pr.PyRanges(
            region_names_to_coordinates(
                region_bin_topics_top_3k[topic].index
            ).sort_values(
                ["Chromosome", "Start", "End"]
            )
        )


dem_region_sets = {
    "top3k": {
        name:
        pr.PyRanges(
            region_names_to_coordinates(
                region_bin_topics_top_3k[name].index
            ).sort_values(
                ["Chromosome", "Start", "End"]
            )
        ) for name in region_bin_topics_top_3k},
    "otsu": {
        name:
        pr.PyRanges(
            region_names_to_coordinates(
                region_bin_topics_otsu[name].index
            ).sort_values(
                ["Chromosome", "Start", "End"]
            )
        ) for name in region_bin_topics_otsu},
    "markers_dict": {
        name:
        pr.PyRanges(
            region_names_to_coordinates(
                markers_dict[name].index
            ).sort_values(
                ["Chromosome", "Start", "End"]
            )
        ) for name in markers_dict}
    }


# Filter multiome data only for the selected regions
df_diff_regions = pd.concat([region_sets[k].dfs[chr]
                             for k in region_sets
                             for chr in region_sets[k].dfs]).drop_duplicates()
df_diff_regions = df_diff_regions.sort_values(by=['Chromosome', 'Start', 'End'])
diff_regions = (df_diff_regions["Chromosome"] +
                ":" + df_diff_regions["Start"].astype(str) +
                "-" + df_diff_regions["End"].astype(str)).values

# Transform as mudata object
# Match RNA barcode
mdata["rna"].obs_names = cell_data.index + "___smpl"
new_mdata = scenicplus.data_wrangling.adata_cistopic_wrangling.process_multiome_data(
    GEX_anndata=mdata["rna"],
    cisTopic_obj=cistopic_obj,
    use_raw_for_GEX_anndata=False,  # ????
    imputed_acc_kwargs=None,
    bc_transform_func=lambda x: x,
    )

pre_atac = new_mdata["scATAC"][:, new_mdata["scATAC"].var.index.isin(
    diff_regions)]
# rename regions
pre_atac.var_names = \
    (pre_atac.var["Chromosome"] +
     "-" + pre_atac.var["Start"].astype(str) +
     "-" + pre_atac.var["End"].astype(str)).values

pre_rna = new_mdata["scRNA"].copy()
pre_atac.layers["counts"] = sp.sparse.csr_matrix(pre_atac.X.copy())
pre_rna.layers["counts"] = sp.sparse.csr_matrix(pre_rna.X.copy())
mdata = mu.MuData({"rna": pre_rna, "atac": pre_atac})
del new_mdata


use_gene_boundaries = args["use_gene_boundaries"]
upstream = tuple([int(num) for num in args["search_space_upstream"].split(' ')])
downstream = tuple([int(num) for num in args["search_space_downstream"].split(' ')])
extend_tss = tuple([int(num) for num in args["search_space_extend_tss"].split(' ')])
remove_promoters = args["remove_promoters"]
importance_scoring_method = args["region_to_gene_importance_method"]
correlation_scoring_method = args["region_to_gene_correlation_method"]
mask_expr_dropout = True

# Download chromosome sizes and gene body coordinates
result = scenicplus.data_wrangling.gene_search_space.download_gene_annotation_and_chromsizes(
        species="hsapiens", #organism,
        biomart_host="http://www.ensembl.org",
        use_ucsc_chromosome_style=True)

if type(result) is tuple:
    annot, chromsizes = result
else:
    annot = result
    print(
        "Chrosomome sizes was not found, please provide this information manually.")

# Calculate search space
mdata["atac"].var_names = mdata["atac"].var_names.str.split('-', 1).str[0]\
    + ':' + mdata["atac"].var_names.str.split('-', 1).str[1]
mdata["rna"][:, mdata["rna"].X.sum(0) != 0]
mdata["atac"][:, mdata["atac"].X.sum(0) != 0]

print(mdata["atac"].var_names)
search_space = scenicplus.data_wrangling.gene_search_space.get_search_space(
        scplus_region=mdata["atac"].var.index,
        scplus_genes=mdata["rna"].var.index,
        gene_annotation=annot,
        chromsizes=chromsizes,
        use_gene_boundaries=use_gene_boundaries,
        upstream=upstream,
        downstream=downstream,
        extend_tss=extend_tss,
        remove_promoters=remove_promoters)

# Calculate regions to genes relationships
p2g = calculate_regions_to_genes_relationships(
        df_exp_mtx=mdata["rna"].to_df(),
        df_acc_mtx=mdata["atac"].to_df(),
        search_space=search_space,
        temp_dir=args["temp_dir"],
        mask_expr_dropout=False,
        importance_scoring_method=importance_scoring_method,
        correlation_scoring_method=correlation_scoring_method,
        n_cpu=njobs,)


dem_max_bg_regions = 500
dem_balance_number_of_promoters = True
dem_promoter_space = 1_000
dem_adj_pval_thr = 0.05
dem_log2fc_thr = 1.0
dem_mean_fg_thr = 0.0
dem_motif_hit_thr = 3.0

# Run DEM
print("Running DEM")
run_motif_enrichment_dem(
    dem_region_sets,
    dem_db_fname=cistarget_score_db_fname,
    output_fname_dem_result=dem_results_path,
    output_fname_dem_html="",
    n_cpu=32,
    path_to_genome_annotation="aertslab/genomes/hg38/hg38_ensdb_v86.csv",
    temp_dir=tmp_scenicplus,
    species=species,
    fraction_overlap_w_dem_database=0.4,
    path_to_motif_annotations=path_to_motif_annotations,
    annotation_version=annotation_version,
    motif_similarity_fdr=0.05,
    orthologous_identity_threshold=0.8,
    annotations_to_use=["Direct_annot", "Orthology_annot"],
    write_html=False,
    max_bg_regions=dem_max_bg_regions,
    balance_number_of_promoters=dem_balance_number_of_promoters,
    promoter_space=dem_promoter_space,
    adjpval_thr=dem_adj_pval_thr,
    log2fc_thr=dem_log2fc_thr,
    mean_fg_thr=dem_mean_fg_thr,
    motif_hit_thr=dem_motif_hit_thr,
    seed=555,
)

# Open cisTarget database
print("Running cistarget")
cistarget_db = pycistarget.motif_enrichment_cistarget.cisTargetDatabase(
    cistarget_ranking_db_fname,
    region_sets=region_sets,
    name="cistarget",
    fraction_overlap=0.4)
# Run cistarget
run_motif_enrichment_cistarget(
    region_sets,
    cistarget_db,
    output_fname_cistarget_result=cistarget_results_path,
    n_cpu=32,
    fraction_overlap_w_cistarget_database=0.4,
    auc_threshold=0.005,
    nes_threshold=3,
    rank_threshold=0.05,
    path_to_motif_annotations=path_to_motif_annotations,
    annotation_version=annotation_version,
    motif_similarity_fdr=0.05,
    orthologous_identity_threshold=0.8,
    temp_dir=tmp_scenicplus,
    species=species,
    annotations_to_use=["Direct_annot", "Orthology_annot"]
)
del cistarget_db

# Reformat cistromes results, merging DEM and Cistarget pairs
# Giving Direct and Extended cistromes
if os.path.exists(dem_results_path):
    paths_to_motif_enrichment_results = [dem_results_path, cistarget_results_path]
else:
    # Save empty hdf file to no thave missiong output for snakemake workflow
    dem_hdf = h5py.File(dem_results_path, 'w')
    dem_hdf.close()

    paths_to_motif_enrichment_results = [cistarget_results_path]
    warnings.warn("Warning...........No significant result from DEM, you might want to adjust thresolds.\nHere, only cistarget method results will be kept.")

output_cistromes_annotations_direct = args["annotation_direct_path"]
output_cistromes_annotations_extended = args["annotation_extended_path"]
output_tf_names = args["tf_names_path"]
prepare_motif_enrichment_results(
    paths_to_motif_enrichment_results=paths_to_motif_enrichment_results,
    multiome_mudata_fname=mdata,
    out_file_direct_annotation=output_cistromes_annotations_direct,
    out_file_extended_annotation=output_cistromes_annotations_extended,
    out_file_tf_names=output_tf_names,
    direct_annotation=["Direct_annot"],
    extended_annotation=["Orthology_annot"])


# Open cistromes direct and extended
direct_h5ad = ad.read_h5ad(output_cistromes_annotations_direct)
direct_h5ad.var["motifs"] = direct_h5ad.var["motifs"].astype(str)+","
extended_h5ad = ad.read_h5ad(output_cistromes_annotations_extended)
extended_h5ad.var["motifs"] = extended_h5ad.var["motifs"].astype(str)+","

# Concat cistromes
## Group motif per TF in var["motifs"]
all_var = pd.concat([direct_h5ad.var, extended_h5ad.var]).reset_index().groupby("index").apply(lambda x: x.sum())
all_var["motifs"] = all_var["motifs"].str.strip(",")
## Sum anndata.X slots
all_h5ad = ad.concat([direct_h5ad, extended_h5ad], join="outer")
all_tfb = all_h5ad.to_df().reset_index().groupby("index").sum()
all_h5ad = ad.AnnData(all_tfb)

with open(output_tf_names) as f:
    tf_names = f.read().split("\n")

# TF to gene relationships
tf_to_gene_prior = infer_TF_to_gene(
    multiome_mudata_fname=mudata_file,
    tf_names=tf_names,
    temp_dir=ray_tmp_dir,
    # adj_out_fname="",#tf_to_gene_prior_path,
    method=args["method_mdl"],
    n_cpu=njobs,
    seed=1)

# Get cistromes
cistromes = all_h5ad

print(p2g.head())
print(tf_to_gene_prior.head())

# Infer eGRN
mdl = infer_grn(
    TF_to_gene_adj=tf_to_gene_prior,
    region_to_gene_adj=p2g,
    cistromes=cistromes,
    ranking_db_fname=ranking_db_fname,
    is_extended=True,
    temp_dir=ray_tmp_dir,
    order_regions_to_genes_by=args.order_regions_to_genes_by,
    order_TFs_to_genes_by=args.order_TFs_to_genes_by,
    gsea_n_perm=args.gsea_n_perm,
    quantiles=args.quantile_thresholds_region_to_gene,
    top_n_regionTogenes_per_gene=args.top_n_regionTogenes_per_gene,
    top_n_regionTogenes_per_region=args.top_n_regionTogenes_per_region,
    binarize_using_basc=True,
    min_regions_per_gene=args.min_regions_per_gene,
    rho_dichotomize_tf2g=True,
    rho_dichotomize_r2g=True,
    rho_dichotomize_eregulon=True,
    keep_only_activating=False,
    rho_threshold=args.rho_threshold,
    min_target_genes=args.min_target_genes,
    n_cpu=1)

mdl = mdl.groupby(["TF", "Gene"])["rho_TF2G"].sum().reset_index()
mdl = mdl.rename({"TF": "source", "Gene": "target", "rho_TF2G": "score"}, axis=1)
mdl.to_csv(output_fname, sep=",", index=False)

