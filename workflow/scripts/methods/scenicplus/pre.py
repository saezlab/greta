import logging
import joblib
import scipy as sp
import os
os.environ['NUMBA_CACHE_DIR'] = '/tmp/'
import pickle
import numpy as np
import scanpy as sc
import muon as mu

# for chromsizes
import pyranges as pr
import requests
import pandas as pd
import pycisTopic
import polars as pl
from pycisTopic.pseudobulk_peak_calling import export_pseudobulk
from typing import Dict
import pathlib
from typing import (List, Literal, TYPE_CHECKING)
import argparse
from pathlib import Path

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--frags', required=True)
parser.add_argument('-i', '--mudata', required=True)
parser.add_argument('-o', '--output', required=True)
parser.add_argument('-t', '--tmp_scenicplus', required=True)
parser.add_argument('-g', '--organism', required=True)
parser.add_argument('-n', '--njobs', required=True, type=int)
parser.add_argument('-m', '--chrom_sizes_m', required=True)
parser.add_argument('-j', '--chrom_sizes_h', required=True)
parser.add_argument('-a', '--annot_human', required=True)
parser.add_argument('-b', '--annot_mouse', required=True)
args = vars(parser.parse_args())

frags = args['frags']
mudata_file = args['mudata']
output_fname = args['output']
tmp_scenicplus = args['tmp_scenicplus']
# ray_tmp_dir = os.path.join(tmp_scenicplus, 'ray_tmp')
ray_tmp_dir = "/local/scratch/tmp"
njobs = args['njobs']
output = args['output']


organism = args['organism']
if organism == 'hg38':
    organism = "hsapiens"
    chromsizes_fname = args['chrom_sizes_h']
    annot_fname = args['annot_human']
elif organism == 'mm10':
    organism = "mmusculus"
    chromsizes_fname = args['chrom_sizes_m']
    annot_fname = args['annot_mouse']

annot = pd.read_csv(annot_fname)
annot.Start = annot.Start.astype(np.int32)
annot.Gene = annot.Gene.astype(str)
print(annot)
annot = pl.DataFrame(annot)
print(annot)

##################
# Set parameters #
##################
dataset_name = 'pbmc10k'

fragments_files = {
    'smpl': frags}

# SCENIC+ paths
# cisTopic paths
cis_topic_tmp_dir = os.path.join(tmp_scenicplus, 'cisTopic')
quality_control_dir = os.path.join(cis_topic_tmp_dir, 'quality_control')
bed_folder = os.path.join(cis_topic_tmp_dir, "pseudobulk_bed_files")
bigwig_folder = os.path.join(cis_topic_tmp_dir, "pseudobulk_bw_files")
bed_pickle = os.path.join(cis_topic_tmp_dir, "pseudobulk_bw_files", 'bed_paths.pkl')
bw_pickle = os.path.join(cis_topic_tmp_dir, "pseudobulk_bw_files", 'bw_paths.pkl')
# MACS paths
macs_folder = os.path.join(cis_topic_tmp_dir, "MACS")
narrow_peaks_pickle = os.path.join(macs_folder, 'narrow_peaks_dict.pkl')
# Consensus peaks paths
consensus_peaks_bed = os.path.join(cis_topic_tmp_dir, 'consensus_region.bed')
quality_control_dir = os.path.join(cis_topic_tmp_dir, 'quality_control')
# cisTopic object path
cistopic_obj_path = os.path.join(cis_topic_tmp_dir, 'cistopic_obj.pkl')

os.makedirs(tmp_scenicplus, exist_ok=True)
os.makedirs(ray_tmp_dir, exist_ok=True)
os.makedirs(cis_topic_tmp_dir, exist_ok=True)
os.makedirs(quality_control_dir, exist_ok=True)
os.makedirs(macs_folder, exist_ok=True)
os.makedirs(bed_folder, exist_ok=True)
os.makedirs(bigwig_folder, exist_ok=True)

# Celltype annotations from MuData object
# Replace spaces with underscores
mudata = mu.read_h5mu(mudata_file)
cell_data = mudata.obs
cell_data['celltype'] = cell_data['celltype'].astype(str) 
cell_data['celltype'] = cell_data['celltype'].str.replace(' ', '_')

# Add cell barcode column
# Create fragments dictionary
fragments_dict = {}
batch_ids = cell_data['batch'].unique()
cell_data['barcode'] = cell_data.index
for batch_id in batch_ids:
    # Add batch_id to fragments_dict
    fragments_dict[batch_id] = fragments_files[batch_id]
    # Remove batch_id from cell barcodes
    cell_data.loc[cell_data['batch'] == batch_id, 'barcode'] = \
        cell_data.loc[cell_data['batch'] == batch_id, 'barcode'].str.replace(
            batch_id + '_', '')


cell_data['barcode'] = cell_data['barcode'].str.split('-1').str[0] + '-1'
cell_data.index = cell_data['barcode']
cell_data.head(3)

# Get chromosome sizes (for hg38 here)
#chromsizes = pd.read_csv(chromsizes_url, sep='\t', header=None)
chromsizes = pd.read_csv(chromsizes_fname, sep='\t', header=None)

chromsizes.columns = ['Chromosome', 'End']
chromsizes['Start'] = [0]*chromsizes.shape[0]
chromsizes = chromsizes.loc[:, ['Chromosome', 'Start', 'End']]
# Exceptionally in this case, to agree with CellRangerARC annotations
chromsizes['Chromosome'] = [chromsizes['Chromosome'][x].replace('v', '.')
                            for x in range(len(chromsizes['Chromosome']))]
chromsizes['Chromosome'] = [chromsizes['Chromosome'][x].split('_')[1]
                            if len(chromsizes['Chromosome'][x].split('_')) > 1
                            else chromsizes['Chromosome'][x]
                            for x in range(len(chromsizes['Chromosome']))]
chromsizes = pr.PyRanges(chromsizes)

bw_paths, bed_paths = export_pseudobulk(
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


# Save all pickle files
pickle.dump(bed_paths,
            open(bed_pickle, 'wb'))
pickle.dump(bw_paths,
            open(bw_pickle, 'wb'))

bed_paths = pickle.load(
    open(os.path.join(bed_pickle), 'rb'))
bw_paths = pickle.load(
    open(os.path.join(bw_pickle), 'rb'))

################
# Peak calling #
################
from pycisTopic.pseudobulk_peak_calling import peak_calling
# Run peak calling
macs_path = 'macs2'
narrow_peaks_dict = peak_calling(macs_path,
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

pickle.dump(narrow_peaks_dict,
            open(narrow_peaks_pickle, 'wb'))

narrow_peaks_dict = pickle.load(
    open(os.path.join(macs_folder, 'narrow_peaks_dict.pkl'), 'rb'))

######################
# Get consensus peaks #
######################
from pycisTopic.iterative_peak_calling import get_consensus_peaks
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
# Enable Polars global string cache so all categoricals are created with the same
# string cache.
pl.enable_string_cache()

from pycisTopic.fragments import (
    get_fragments_in_peaks,
    get_fragments_per_cb,
    get_insert_size_distribution,
)
from pycisTopic.topic_binarization import threshold_otsu
from pycisTopic.tss_profile import get_tss_profile
from scipy.stats import gaussian_kde

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


from pycisTopic.qc import get_barcodes_passing_qc_for_sample
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
    #write barcodes passing QC to file
    pd.Series(sample_id_to_barcodes_passing_filters[sample_id]).to_csv(
        os.path.join(quality_control_dir, f'{sample_id}.barcodes_passing_qc.txt'),
        header=None, index=None, sep='\t'
    )

fragments_df_pl = pycisTopic.fragments.read_bed_to_polars_df(
        bed_filename=fragments_dict[sample_id],
        engine="polars",
        min_column_count=4,
    ).lazy()

from pycisTopic.cistopic_class import create_cistopic_object_from_fragments
import polars as pl

cistopic_obj_list = []
for sample_id in fragments_dict:
    # compute matrix from fragments, peaks, and barcodes passing QC.
    counts_fragments_matrix, cbs, region_ids = create_fragment_matrix_from_fragments(
        fragments_bed_filename=fragments_dict[sample_id],
        regions_bed_filename=consensus_peaks_bed,
        barcodes_tsv_filename=os.path.join(quality_control_dir, f'{sample_id}.barcodes_passing_qc.txt'),
        # blacklist_bed_filename="hg38-blacklist.v2.bed",  # Or None
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

pickle.dump(
    cistopic_obj,
    open(os.path.join(cis_topic_tmp_dir, "cistopic_obj.pkl"), "wb")
)

cistopic_obj = pickle.load(
    open(os.path.join(cis_topic_tmp_dir, "cistopic_obj.pkl"), "rb")
)

import pycisTopic.lda_models
# Run models
models = pycisTopic.lda_models.run_cgs_models(
    cistopic_obj,
    n_topics=[2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50],
    n_cpu=njobs,
    n_iter=150,
    random_state=555,
    alpha=50,
    alpha_by_topic=True,
    eta=0.1,
    eta_by_topic=False,
    save_path=cis_topic_tmp_dir,
    _temp_dir=ray_tmp_dir
)

pickle.dump(
    models,
    open(os.path.join(cis_topic_tmp_dir, "models.pkl"), "wb")
)

model = pycisTopic.lda_models.evaluate_models(
    models,
    select_model=None,
    return_model=True
)
cistopic_obj.add_LDA_model(model)
pickle.dump(
    cistopic_obj,
    open(os.path.join(
        cis_topic_tmp_dir, "cistopic_obj.pkl"), "wb")
)

pickle.load(
    open(os.path.join(
        cis_topic_tmp_dir, "cistopic_obj.pkl"), "rb")
)

# Now, we keep only variable regions
from pycisTopic.diff_features import (
    impute_accessibility,
    normalize_scores,
    find_highly_variable_features,
    find_diff_features
)
from pycisTopic.topic_binarization import binarize_topics

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

from pycisTopic.utils import region_names_to_coordinates
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

# Transform as mudata object
# Match RNA barcode
mudata["rna"].obs_names = cell_data.index + "___smpl"
import scenicplus.data_wrangling.adata_cistopic_wrangling
new_mudata = scenicplus.data_wrangling.adata_cistopic_wrangling.process_multiome_data(
    GEX_anndata=mudata["rna"],
    cisTopic_obj=cistopic_obj,
    use_raw_for_GEX_anndata=False, # ????
    imputed_acc_kwargs=None,
    bc_transform_func=lambda x: x,
    )

pre_atac = new_mudata["scATAC"][:, new_mudata["scATAC"].var.index.isin(
    diff_regions)]
# rename regions
pre_atac.var_names = \
    (pre_atac.var["Chromosome"] +
     "-" + pre_atac.var["Start"].astype(str) +
     "-" + pre_atac.var["End"].astype(str)).values
pre_mudata = mu.MuData({"rna": new_mudata["scRNA"], "atac": pre_atac})

# Add counts layers
pre_mudata["atac"].layers["counts"] = pre_mudata["atac"].X.copy()
pre_mudata["rna"].layers["counts"] = mudata["rna"].X.copy()

pre_mudata.write_h5mu(output_fname)
