# Rerun cistopic, needed for cistarget
import matplotlib.pyplot as plt
from typing import List, Literal, Optional, Tuple, Set, Dict, Union
import argparse
import pathlib
import numpy as np
import os
import joblib
os.environ["NUMBA_CACHE_DIR"] = "/tmp"
import muon as mu
import pandas as pd
import polars as pl
import pyranges as pr
import pycisTopic.cistopic_class
from pycisTopic.topic_binarization import binarize_topics
from pycisTopic.utils import region_names_to_coordinates
from pycisTopic.diff_features import (
    impute_accessibility,
    normalize_scores,
    find_highly_variable_features,
    find_diff_features
)
from scenicplus.cli.commands import (
    prepare_motif_enrichment_results,
    _run_dem_single_region_set,
    _get_foreground_background
)
import pycistarget.motif_enrichment_cistarget


parser = argparse.ArgumentParser()
parser.add_argument('-i', '--mudata', required=True)
parser.add_argument('-p', '--p2g', required=True)
parser.add_argument('-d', '--raw_data', required=True)
parser.add_argument('-t', '--cistarget_results', required=True)
parser.add_argument('-u', '--dem_results', required=True)
parser.add_argument('-o', '--output', required=True)
parser.add_argument('-g', '--organism', required=True)
parser.add_argument('-s', '--dem_scores', required=True)
parser.add_argument('-r', '--cistarget_rankings', required=True)
parser.add_argument('-c', '--njobs', required=True, type=int)
parser.add_argument('-a', '--annotation_direct_path', required=True)
parser.add_argument('-b', '--annotation_extended_path', required=True)
parser.add_argument('-n', '--tf_names_path', required=True)
args = parser.parse_args()

mudata_path = args.mudata
cistarget_rankings_fname = args.cistarget_rankings
dem_scores_fname = args.dem_scores
cistarget_results_path = args.cistarget_results
dem_results_path = args.dem_results
raw_data = args.raw_data
p2g = args.p2g
organism = args.organism
njobs = args.njobs

output_cistromes_annotations_direct = args.annotation_direct_path
output_cistromes_annotations_extended = args.annotation_extended_path
output_tf_names = args.tf_names_path

if organism == "hg38":
    annotation_version = "hg38"
elif organism == "mm10":
    annotation_version = "mm10"
else:
    raise ValueError("Invalid organism")

# Keep peaks contained in enhancers only
p2g = pd.read_csv(p2g)
regions = p2g['cre'].unique()

# Load the (raw) data and filter regions
mudata = mu.read_h5mu(mudata_path)
# mask columns in regions
mask = mudata["atac"].var_names.isin(regions)
cistopic_obj = pycisTopic.cistopic_class.create_cistopic_object(
    mudata["atac"].layers["counts"][:, mask].T.toarray(),
    cell_names=mudata["atac"][:, mask].obs_names.values.tolist(),
    region_names=mudata["atac"][:, mask].var_names.values.tolist()
    )

# Load the celltype annotations
mudata = mu.read_h5mu(raw_data)
cell_data = mudata.obs
cell_data['celltype'] = cell_data['celltype'].astype(str) 
cell_data['celltype'] = cell_data['celltype'].str.replace(' ', '_')
cistopic_obj.add_cell_data(cell_data, split_pattern='-')

# Rerun cistopic
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
    save_path="/tmp",
    _temp_dir="/tmp"
)
model = pycisTopic.lda_models.evaluate_models(
    models,
    select_model=None,
    return_model=True
)
cistopic_obj.add_LDA_model(model)

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
    _temp_dir="/tmp",
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
    "top3k": region_bin_topics_top_3k,
    "otsu": region_bin_topics_otsu,
    "markers_dict": markers_dict
    }

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
        temp_folder="/tmp"
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
    region_set_dict: Dict[str, pr.PyRanges] = {}

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
        for name, foreground_region_sets, background_region_sets
        in _get_foreground_background(region_set_dict)
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
            dem_result.write_hdf5(
                path = output_fname_dem_result,
                mode = "a"
            )


dem_max_bg_regions = 500
dem_balance_number_of_promoters = True
dem_promoter_space = 1_000
dem_adj_pval_thr = 0.05
dem_log2fc_thr = 1.0
dem_mean_fg_thr = 0.0
dem_motif_hit_thr = 3.0
# Run DEM
run_motif_enrichment_dem(
    region_sets,
    dem_db_fname=dem_scores_fname,
    output_fname_dem_result=dem_results_path,
    output_fname_dem_html="",
    n_cpu=1,
    temp_dir="/tmp",
    species=organism,
    fraction_overlap_w_dem_database=0.4,
    path_to_motif_annotations="path_to_motif_annotations",
    annotation_version=annotation_version,
    motif_similarity_fdr=0.05,
    orthologous_identity_threshold=0.8,
    annotations_to_use=["motif", "orthologous"],
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

# Run cisTarget
cistarget_db = pycistarget.motif_enrichment_cistarget.cisTargetDatabase(
    cistarget_rankings_fname,
    region_sets=region_sets,
    name="cistarget",
    fraction_overlap=0.4)

# Run cistarget
run_motif_enrichment_cistarget(
    region_sets,
    cistarget_db,
    save_path=os.path.join(cistarget_results_path),
    n_cpu=1,
    fraction_overlap_w_cistarget_database=0.4,
    auc_threshold=0.005,
    nes_threshold=3,
    rank_threshold=0.05,
    path_to_motif_annotations="path_to_motif_annotations",
    annotation_version=annotation_version,
    motif_similarity_fdr=0.05,
    orthologous_identity_threshold=0.8,
    temp_dir="/tmp",
    species=organism,
    annotations_to_use=["motif", "orthologous"]
)

prepare_motif_enrichment_results(
    paths_to_motif_enrichment_results=[dem_results_path, cistarget_results_path],
    multiome_mudata_fname=mudata_path,
    out_file_direct_annotation=output_cistromes_annotations_direct,
    out_file_extended_annotation=output_cistromes_annotions_extended,
    out_file_tf_names=output_tf_names,
    direct_annotation="Direct_annot",
    extended_annotation="Orthology_annot")
