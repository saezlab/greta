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

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--mudata', required=True)
parser.add_argument('-p', '--p2g', required=True)
parser.add_argument('-o', '--output', required=True)
parser.add_argument('-g', '--organism', required=True)
parser.add_argument('-s', '--cistarget_scores', required=True)
parser.add_argument('-r', '--cistarget_rankings', required=True)
parser.add_argument('-c', '--njobs', required=True, type=int)
args = parser.parse_args()

mudata_path = args.mudata
cistarget_scores_fname = args.cistarget_scores
cistarget_ranking_fname = args.cistarget_rankings
cistarget_results_path = args.output
organism = args.organism
njobs = args.njobs

if organism == "human":
    annotation_version = "hg38"
elif organism == "mouse":
    annotation_version = "mm10"

# Load the data
mudata = mu.read_h5mu(mudata_path)
cistopic_obj = pycisTopic.cistopic_class.create_cistopic_object(
    mudata["atac"].X.T,
    cell_names=mudata["atac"].obs_names.values.tolist(),
    region_names=mudata["atac"].var_names.values.tolist()
    )

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

# Run cisTarget
import pycistarget.motif_enrichment_cistarget
cistarget_db = pycistarget.motif_enrichment_cistarget.cisTargetDatabase(
    cistarget_ranking_fname,
    region_sets=region_sets,
    name="cistarget",
    fraction_overlap=0.4)


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


run_motif_enrichment_cistarget(
    region_sets,
    cistarget_db,
    os.path.join(cistarget_results_path),
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


#pd.read_hdf(os.path.join(scenicplus_tmp_dir, "cistarget_results.h5"), key='cistarget_results')
