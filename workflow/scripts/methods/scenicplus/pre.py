import scipy as sp
import os
os.environ['NUMBA_CACHE_DIR'] = '/tmp/'
import pickle
import numpy as np
import scanpy as sc
import muon as mu
import pyranges as pr
import pandas as pd
import pycisTopic
import polars as pl
from pycisTopic.pseudobulk_peak_calling import export_pseudobulk
import argparse
from pathlib import Path
from pycisTopic.pseudobulk_peak_calling import peak_calling
from pycisTopic.iterative_peak_calling import get_consensus_peaks
from pycisTopic.cistopic_class import create_cistopic_object_from_fragments
from pycisTopic.diff_features import (
    impute_accessibility,
    normalize_scores,
    find_highly_variable_features,
    find_diff_features
)
import h5py
import pycisTopic.lda_models
from pycisTopic.topic_binarization import binarize_topics
from pycisTopic.utils import region_names_to_coordinates
import scenicplus.data_wrangling.adata_cistopic_wrangling
import importlib
import logging
from pycisTopic.fragments import (
    get_fragments_in_peaks,
    get_fragments_per_cb,
    get_insert_size_distribution
)
from pycisTopic.topic_binarization import threshold_otsu
from pycisTopic.tss_profile import get_tss_profile
from scipy.stats import gaussian_kde


def compute_qc_stats(
    fragments_df_pl,
    regions_df_pl,
    tss_annotation,
    tss_flank_window=2000,
    tss_smoothing_rolling_window=10,
    tss_minimum_signal_window=100,
    tss_window=50,
    tss_min_norm=0.2,
    use_genomic_ranges=True,
    min_fragments_per_cb=10,
    collapse_duplicates=True,
    no_threads=8,
):

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
    logger.info("Get basic fragments statistics per cell barcode.")
    fragments_stats_per_cb_df_pl = get_fragments_per_cb(
        fragments_df_pl=fragments_df_pl,
        min_fragments_per_cb=min_fragments_per_cb,
        collapse_duplicates=collapse_duplicates,
    ).rename({"by":"CB"})
    logger.info("Get total fragment counts and unique fragment counts per region.")
    fragments_in_peaks_df_pl = get_fragments_in_peaks(
        fragments_df_pl=fragments_df_pl,
        regions_df_pl=regions_df_pl,
    ).rename({"by":"CB"})
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
    ).collect()
    logger.info("Get insert size distribution of fragments.")
    insert_size_dist_df_pl = get_insert_size_distribution(
        fragments_df_pl=fragments_df_pl,
    )
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
    logger.info("Add TSS enrichment to fragments statistics per cell barcode.")
    fragments_stats_per_cb_df_pl = (
        fragments_stats_per_cb_df_pl.lazy().join(
            tss_enrichment_per_cb.lazy(),
            how="left",
            on="CB",
        )
        .with_columns(
            pl.col("tss_enrichment").fill_null(0.0),
        )
        .collect()
    )
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
    kde_data_for_tss_enrichment = np.vstack(
        [log10_unique_fragments_in_peaks_count, tss_enrichment]
    )
    kde_data_for_fraction_of_fragments_in_peaks = np.vstack(
        [log10_unique_fragments_in_peaks_count, fraction_of_fragments_in_peaks]
    )
    kde_data_for_duplication_ratio = np.vstack(
        [log10_unique_fragments_in_peaks_count, duplication_ratio]
    )
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
    fragments_tsv_filename,
    regions_bed_filename,
    tss_annotation_bed_filename,
    output_prefix,
    tss_flank_window=2000,
    tss_smoothing_rolling_window=10,
    tss_minimum_signal_window=100,
    tss_window=50,
    tss_min_norm=0.2,
    use_genomic_ranges=True,
    min_fragments_per_cb=10,
    collapse_duplicates=True,
    no_threads=8,
):
    import logging

    from pycisTopic.fragments import read_bed_to_polars_df, read_fragments_to_polars_df
    from pycisTopic.gene_annotation import read_tss_annotation_from_bed
    from pycisTopic.qc import compute_qc_stats, get_otsu_threshold

    # Remove trailing dot(s) from the output prefix.
    output_prefix = output_prefix.rstrip(".")

    class RelativeSeconds(logging.Formatter):
        def format(self, record):
            record.relativeCreated = record.relativeCreated // 1000
            return super().format(record)

    formatter = RelativeSeconds(
        "%(asctime)s %(relativeCreated)ds - %(levelname)s - %(name)s:%(funcName)s - %(message)s"
    )
    logging.basicConfig(
        # format=formatter,
        filename=f"{output_prefix}.pycistopic_qc.log",
        encoding="utf-8",
        level=logging.INFO,
    )
    logging.root.handlers[0].setFormatter(formatter)

    logger = logging.getLogger(__name__)

    logger.info(f'Loading TSS annotation from "{tss_annotation_bed_filename}".')
    tss_annotation_bed_df_pl = read_tss_annotation_from_bed(
        tss_annotation_bed_filename=tss_annotation_bed_filename
    )

    logger.info(f'Loading regions BED file from "{regions_bed_filename}".')
    regions_df_pl = read_bed_to_polars_df(
        bed_filename=regions_bed_filename,
        min_column_count=3,
    )

    logger.info(f'Loading fragments TSV file from "{fragments_tsv_filename}".')
    fragments_df_pl = read_fragments_to_polars_df(
        fragments_tsv_filename,
        engine="pyarrow",
    )

    logger.info("Computing QC stats.")
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

    logger.info(f'Writing "{output_prefix}.fragments_stats_per_cb.parquet".')
    fragments_stats_per_cb_df_pl.write_parquet(
        f"{output_prefix}.fragments_stats_per_cb.parquet",
        compression="zstd",
        use_pyarrow=True,
    )

    logger.info(f'Writing "{output_prefix}.fragments_insert_size_dist.parquet".')
    insert_size_dist_df_pl.write_parquet(
        f"{output_prefix}.fragments_insert_size_dist.parquet",
        compression="zstd",
        use_pyarrow=True,
    )

    logger.info(f'Writing "{output_prefix}.tss_norm_matrix_sample.parquet".')
    tss_norm_matrix_sample.write_parquet(
        f"{output_prefix}.tss_norm_matrix_sample.parquet",
        compression="zstd",
        use_pyarrow=True,
    )

    logger.info(f'Writing "{output_prefix}.tss_norm_matrix_per_cb.parquet".')
    tss_norm_matrix_per_cb.write_parquet(
        f"{output_prefix}.tss_norm_matrix_per_cb.parquet",
        compression="zstd",
        use_pyarrow=True,
    )

    logger.info("Calculating Otsu thresholds.")
    (
        unique_fragments_in_peaks_count_otsu_threshold,
        tss_enrichment_otsu_threshold,
        fragments_stats_per_cb_for_otsu_threshold_df_pl,
    ) = get_otsu_threshold(
        fragments_stats_per_cb_df_pl=fragments_stats_per_cb_df_pl,
        min_otsu_fragments=100,
        min_otsu_tss=1.0,
    )

    logger.info(
        f'Writing "{output_prefix}.fragments_stats_per_cb_for_otsu_thresholds.parquet".'
    )
    fragments_stats_per_cb_for_otsu_threshold_df_pl.write_parquet(
        f"{output_prefix}.fragments_stats_per_cb_for_otsu_thresholds.parquet",
        compression="zstd",
        use_pyarrow=True,
    )
    logger.info(
        f'Writing "{output_prefix}.fragments_stats_per_cb_for_otsu_thresholds.tsv".'
    )
    fragments_stats_per_cb_for_otsu_threshold_df_pl.write_csv(
        f"{output_prefix}.fragments_stats_per_cb_for_otsu_thresholds.tsv",
        separator="\t",
        include_header=True,
    )

    logger.info(f'Writing "{output_prefix}.cbs_for_otsu_thresholds.tsv".')
    fragments_stats_per_cb_for_otsu_threshold_df_pl.select(pl.col("CB")).write_csv(
        f"{output_prefix}.cbs_for_otsu_thresholds.tsv",
        separator="\t",
        include_header=False,
    )

    logger.info(f'Writing "{output_prefix}.otsu_thresholds.tsv".')
    with open(f"{output_prefix}.otsu_thresholds.tsv", "w") as fh:
        print(
            "unique_fragments_in_peaks_count_otsu_threshold\ttss_enrichment_otsu_threshold\n"
            f"{unique_fragments_in_peaks_count_otsu_threshold}\t{tss_enrichment_otsu_threshold}",
            file=fh,
        )
    logger.info("pycisTopic QC finished.")

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
    fragments_bed_filename,
    engine='polars',
    sample_id=None,
    cb_end_to_remove="-1",
    cb_sample_separator="___",
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


pycisTopic.qc = importlib.import_module("pycisTopic.qc")
pycisTopic.QC = importlib.import_module("pycisTopic.cli.subcommand.qc")
pycisTopic.qc.compute_qc_stats = compute_qc_stats
pycisTopic.QC.qc = qc
get_fragments_per_cb = pycisTopic.qc.get_fragments_per_cb
get_fragments_in_peaks = pycisTopic.qc.get_fragments_in_peaks
compute_kde = pycisTopic.qc.compute_kde


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--frags', nargs='+', required=True)
parser.add_argument('-i', '--mudata', required=True)
parser.add_argument('-o', '--output', required=True)
parser.add_argument('-t', '--tmp_scenicplus', required=True)
parser.add_argument('-z', '--ray_tmp_dir', required=True)
parser.add_argument('-g', '--organism', required=True)
parser.add_argument('-n', '--njobs', required=True, type=int)
parser.add_argument('-a', '--gannot', required=True, nargs='+')
parser.add_argument('-b', '--csize', required=True, nargs='+')
parser.add_argument('-c', '--cistopic_obj', required=True)
args = vars(parser.parse_args())

frags = args['frags']
mudata_file = args['mudata']
output_fname = args['output']
tmp_scenicplus = args['tmp_scenicplus']
ray_tmp_dir = os.path.join(args['ray_tmp_dir'])
njobs = args['njobs']
output = args['output']
gannot = np.array(args['gannot'])
csize = np.array(args['csize'])

# Create tmp dirs
cis_topic_tmp_dir = os.path.join(tmp_scenicplus, 'cisTopic')
quality_control_dir = os.path.join(cis_topic_tmp_dir, 'quality_control')
bed_folder = os.path.join(cis_topic_tmp_dir, "pseudobulk_bed_files")
bigwig_folder = os.path.join(cis_topic_tmp_dir, "pseudobulk_bw_files")
macs_folder = os.path.join(cis_topic_tmp_dir, "MACS")
consensus_peaks_bed = os.path.join(cis_topic_tmp_dir, 'consensus_region.bed')
cistopic_obj_path = os.path.join(args['cistopic_obj'])
os.makedirs(tmp_scenicplus, exist_ok=True)
os.makedirs(quality_control_dir, exist_ok=True)
os.makedirs(cis_topic_tmp_dir, exist_ok=True)
os.makedirs(bed_folder, exist_ok=True)
os.makedirs(bigwig_folder, exist_ok=True)
os.makedirs(macs_folder, exist_ok=True)

# Extract chrom sizes and gannot for given organism
organism = args['organism']
msk = np.array([os.path.basename(p).replace('cist_', '').startswith(organism) for p in gannot])
chromsizes_fname = csize[msk][0]
annot_fname = gannot[msk][0]

annot = pd.read_csv(annot_fname, sep="\t")
annot.Start = annot.Start.astype(np.int32)
annot.Gene = annot.Gene.astype(str)
annot = pl.DataFrame(annot)

# Read metadata
def extract_cat(f, col):
    codes = f[col]['codes'][:]
    cats = f[col]['categories'][:].astype('U')
    return cats[codes]
with h5py.File(mudata_file, "r") as f:
    cell_data = pd.DataFrame(index=f['obs']['_index'][:].astype('U'))
    cell_data['celltype'] = extract_cat(f['obs'], 'celltype')
    cell_data['batch'] = extract_cat(f['obs'], 'batch')
cell_data['celltype'] = cell_data['celltype'].str.replace(' ', '_')
cell_data['barcode'] = cell_data.index

# Make scATAC psbulks
chromsizes = pr.PyRanges(pd.read_csv(chromsizes_fname, sep='\t'))
fragments_dict = {os.path.basename(frag).replace('.frags.tsv.gz', ''): frag for frag in frags}
bw_paths, bed_paths = export_pseudobulk(
    input_data=cell_data,
    variable='celltype',
    sample_id_col='batch',
    chromsizes=chromsizes,
    bed_path=bed_folder,
    bigwig_path=bigwig_folder,
    path_to_fragments=fragments_dict,
    n_cpu=4,
    normalize_bigwig=True,
    temp_dir=ray_tmp_dir,
    split_pattern='_'
)

# Run peak calling
if organism == 'hg38':
    genome_size = 'hs'
elif organism == 'mm10':
    genome_size = 'mm'

narrow_peaks_dict = peak_calling(
    'macs2',
    bed_paths,
    macs_folder,
    input_format='BED',
    genome_size=genome_size,
    n_cpu=njobs,
    _temp_dir=ray_tmp_dir
)


# Get consensus peaks
consensus_peaks = get_consensus_peaks(
    narrow_peaks_dict,
    peak_half_width=250,
    chromsizes=chromsizes,
    path_to_blacklist=None
)
consensus_peaks.to_bed(
    path=consensus_peaks_bed,
    keep=True,
    compression='infer',
    chain=False
)

# Run QC
for key in fragments_dict.keys():
    pycisTopic.QC.qc(
        fragments_tsv_filename=fragments_dict[key],
        regions_bed_filename=consensus_peaks_bed,
        tss_annotation_bed_filename=annot_fname,
        output_prefix=os.path.join(quality_control_dir, key),
        no_threads=njobs,
    )
    
# Get barcodes passing QC
sample_id_to_barcodes_passing_filters = {}
sample_id_to_thresholds = {}
for sample_id in fragments_dict:
    (
        sample_id_to_barcodes_passing_filters[sample_id],
        sample_id_to_thresholds[sample_id]
    ) = pycisTopic.qc.get_barcodes_passing_qc_for_sample(
            sample_id=sample_id,
            pycistopic_qc_output_dir=quality_control_dir,
            use_automatic_thresholds=True,
    )

cistopic_obj_list = []
for sample_id in fragments_dict:
    # compute matrix from fragments, peaks, and barcodes passing QC.
    counts_fragments_matrix, cbs, region_ids = create_fragment_matrix_from_fragments(
        fragments_bed_filename=fragments_dict[sample_id],
        regions_bed_filename=consensus_peaks_bed,
        barcodes_tsv_filename=os.path.join(quality_control_dir, f'{sample_id}.cbs_for_otsu_thresholds.tsv'),
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
cistopic_obj.add_cell_data(cell_data)
cistopic_obj.cell_data.index = cistopic_obj.cell_data.index.str.split('___', expand=True).get_level_values(0)
cistopic_obj.cell_names = cistopic_obj.cell_data.index
fbarcodes = cistopic_obj.cell_data.dropna().index
cistopic_obj.subset(cells=fbarcodes)

pickle.dump(
    cistopic_obj,
    open(os.path.join(quality_control_dir, "cistopic_obj.pkl"), "wb")
)

# Run models
n_topics=[2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
models = pycisTopic.lda_models.run_cgs_models(
    cistopic_obj,
    n_topics=n_topics,
    n_cpu=njobs,
    save_path=cis_topic_tmp_dir,
    _temp_dir=ray_tmp_dir,
)

model = pycisTopic.lda_models.evaluate_models(
    models,
    select_model=None,
    return_model=True
)
cistopic_obj.add_LDA_model(model)
pickle.dump(
    cistopic_obj,
    open(cistopic_obj_path, "wb")
)

obsm = dict()
with h5py.File(mudata_file, "r") as f:
    obsm['X_spectral'] = f['obsm']['X_spectral'][:]
    obsm['X_umap'] = f['obsm']['X_umap'][:]

for k in obsm:
    dim = pd.DataFrame(obsm[k], index=cell_data.index, columns=[f'{k}__{i}' for i in range(obsm[k].shape[1])])
    cisTopic_obj.projections["cell"][k] = dim.loc[fbarcodes, :]

# Find variable regions
imputed_acc_obj = impute_accessibility(
    cistopic_obj,
    selected_cells=None,
    selected_regions=None,
    scale_factor=10**6
)
normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)
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
    split_pattern='_'
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

# Filter multiome data only for the selected regions
df_diff_regions = pd.concat([region_sets[k].dfs[chr]
                             for k in region_sets
                             for chr in region_sets[k].dfs]).drop_duplicates()
df_diff_regions = df_diff_regions.sort_values(by=['Chromosome', 'Start', 'End'])
diff_regions = (df_diff_regions["Chromosome"] +
                ":" + df_diff_regions["Start"].astype(str) +
                "-" + df_diff_regions["End"].astype(str)).values

new_mudata = scenicplus.data_wrangling.adata_cistopic_wrangling.process_multiome_data(
    GEX_anndata=mu.read(os.path.join(mudata_file, 'rna')),
    cisTopic_obj=cistopic_obj,
    use_raw_for_GEX_anndata=False,
    imputed_acc_kwargs=None,
    bc_transform_func=lambda x: x,
    )

print(new_mudata.obs_names)
pre_atac = new_mudata["scATAC"][:, new_mudata["scATAC"].var.index.isin(
    diff_regions)]
# rename regions
pre_atac.var_names = \
    (pre_atac.var["Chromosome"] +
     "-" + pre_atac.var["Start"].astype(str) +
     "-" + pre_atac.var["End"].astype(str)).values

pre_rna = new_mudata["scRNA"].copy()
pre_atac.layers["counts"] = sp.sparse.csr_matrix(pre_atac.X.copy())
pre_rna.layers["counts"] = sp.sparse.csr_matrix(pre_rna.X.copy())


pre_mudata = mu.MuData({"rna": pre_rna, "atac": pre_atac})


pre_mudata.write_h5mu(output_fname)
