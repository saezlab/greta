import joblib
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
from typing import List
from typing import Literal
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--frags', required=True)
parser.add_argument('-m', '--mudata', required=True)
parser.add_argument('-d', '--output', required=True)
parser.add_argument('-t', '--tmp_scenicplus', required=True)
parser.add_argument('-g', '--organism', required=True)
parser.add_argument('-n', '--njobs', required=True)
args = vars(parser.parse_args())

frags = args['frags']
mudata = args['mudata']
output = args['output']
tmp_scenicplus = args['tmp_scenicplus']
ray_tmp_dir = os.path.join(tmp_scenicplus, 'ray_tmp')
organism = args['organism']
njobs = args['njobs']

##################
# Set parameters #
##################
organism = 'hsapiens'
dataset_name = 'pbmc10k'

# Set file names
mudata_file = os.path.join(
    "datasets/{dataset}/annotated.h5mu".format(
        dataset=dataset_name))
fragments_files = {
    'smpl':
    os.path.join("datasets/{dataset}_smpl_fragments.tsv.gz".format(
        dataset=dataset_name))}
# Set temp output file
output_filename = os.path.join(tmp_scenicplus, 'scenicplus.h5mu')

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

if organism == 'hsapiens':
    chromsizes_url = 'http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes'
if organism == 'mmusculus':
    chromsizes_url = 'http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes'


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
chromsizes = pd.read_csv(chromsizes_url, sep='\t', header=None)
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
import pybiomart as pbm
if organism == 'hsapiens':
    dataset = pbm.Dataset(name='hsapiens_gene_ensembl',  host='http://www.ensembl.org')

if organism == 'mmusculus':
    dataset = pbm.Dataset(name='mmusculus_gene_ensembl',  host='http://www.ensembl.org')

annot = dataset.query(
    attributes=[
        'chromosome_name',
        'transcription_start_site',
        'strand', 'external_gene_name',
        'transcript_biotype'
        ])
annot['Chromosome/scaffold name'] = annot['Chromosome/scaffold name'].to_numpy(dtype = str)
filter = annot['Chromosome/scaffold name'].str.contains('CHR|GL|JH|MT')
annot = annot[~filter]
annot['Chromosome/scaffold name'] = annot['Chromosome/scaffold name'].str.replace(r'(\b\S)', r'chr\1')
annot.columns=['Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
annot = annot[annot.Transcript_type == 'protein_coding']
annot["Strand"] = annot["Strand"].replace({1: "+", -1: "-"})
annot.Start = annot.Start.astype(np.int32)
annot = pl.DataFrame(annot)


from __future__ import annotations
import os
from typing import TYPE_CHECKING, Literal
import polars as pl
from pathlib import Path
# Enable Polars global string cache so all categoricals are created with the same
# string cache.
pl.enable_string_cache()


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
    from pycisTopic.qc import compute_qc_stats, get_otsu_threshold
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
    ) = pycisTopic.qc.compute_qc_stats(
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

from pycisTopic.cistopic_class import create_cistopic_object_from_fragments
import polars as pl

cistopic_obj_list = []
for sample_id in fragments_dict:
    sample_metrics = pl.read_parquet(
        os.path.join(
            quality_control_dir, f'{sample_id}.fragments_stats_per_cb.parquet')
    ).to_pandas().set_index("CB").loc[sample_id_to_barcodes_passing_filters[sample_id]]
    cistopic_obj = create_cistopic_object_from_fragments(
        path_to_fragments=fragments_dict[sample_id],
        path_to_regions=consensus_peaks_bed,
        # = path_to_blacklist,
        metrics=sample_metrics,
        valid_bc=sample_id_to_barcodes_passing_filters[sample_id],
        n_cpu=1,
        project=sample_id,
        split_pattern='-'
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
    save_path=tmp_scenicplus,
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
markers_dict= find_diff_features(
    cistopic_obj,
    imputed_acc_obj,
    variable='celltype',
    var_features=variable_regions,
    contrasts=None,
    adjpval_thr=0.05,
    log2fc_thr=np.log2(1.5),
    n_cpu=5,
    _temp_dir=ray_tmp_dir,
    split_pattern = '-'
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
mudata["rna"].obs_names = cell_data.index + "-smpl___smpl"
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

pre_mudata.write_h5mu(output_filename)
