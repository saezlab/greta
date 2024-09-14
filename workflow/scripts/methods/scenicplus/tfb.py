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
import mudata
from scenicplus.cli.commands import (
#    prepare_motif_enrichment_results,
    _get_foreground_background
)
import pycistarget.motif_enrichment_cistarget

import logging
log = logging.getLogger("SCENIC+")

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
parser.add_argument('-z', '--path_to_motif_annotations_mouse')
parser.add_argument('-y', '--path_to_motif_annotations_human')
parser.add_argument('-k', '--temp_dir')
parser.add_argument('-j', '--ray_tmp_dir')
args = parser.parse_args()

mudata_path = args.mudata
cistarget_rankings_fname = args.cistarget_rankings
dem_scores_fname = args.dem_scores
cistarget_results_path = args.cistarget_results
dem_results_path = args.dem_results
raw_data = args.raw_data
p2g = args.p2g
output_tfb = args.output
organism = args.organism
njobs = args.njobs
temp_dir = args.temp_dir
ray_tmp_dir = args.ray_tmp_dir

output_cistromes_annotations_direct = args.annotation_direct_path
output_cistromes_annotations_extended = args.annotation_extended_path
output_tf_names = args.tf_names_path

if organism == "hg38":
    species = "homo_sapiens",
    annotation_version = "v10nr_clust"
    path_to_motif_annotations = args.path_to_motif_annotations_human
elif organism == "mm10":
    species = "mus_musculus"
    annotation_version = "v10nr_clust"
    path_to_motif_annotations = args.path_to_motif_annotations_mouse
else:
    raise ValueError("Invalid organism")

# Keep peaks contained in enhancers only
p2g = pd.read_csv(p2g)
regions = p2g['cre'].unique()

# Load the (raw) data and filter regions
mudata_file = mu.read_h5mu(mudata_path)
# mask columns in regions
mask = mudata_file["atac"].var_names.isin(regions)
cistopic_obj = pycisTopic.cistopic_class.create_cistopic_object(
    mudata_file["atac"].layers["counts"][:, mask].T.toarray(),
    cell_names=mudata_file["atac"][:, mask].obs_names.values.tolist(),
    region_names=mudata_file["atac"][:, mask].var_names.values.tolist()
    )

# Load the celltype annotations
rawdata_file = mu.read_h5mu(raw_data)
cell_data = rawdata_file.obs
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
    save_path=temp_dir,
    _temp_dir=ray_tmp_dir
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
    _temp_dir=ray_tmp_dir,
    split_pattern='-'
)

region_bin_topics_otsu = binarize_topics(
    cistopic_obj, method='otsu',
    plot=True, num_columns=5
)

n_top = int(min(3_000, len(mudata_file["atac"].var_names)/5))
region_bin_topics_top_3k = binarize_topics(
    cistopic_obj, method='ntop', ntop=n_top,
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
            dem_result.write_hdf5(
                path = output_fname_dem_result,
                mode = "a"
            )
from pycistarget.motif_enrichment_dem import (
    # DEM,
    ranksums_numba_multiple,
    mean_axis1, get_log2_fc,
    p_adjust_bh, get_optimal_threshold_roc
    )

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
        orthologous_identity_threshold) -> DEM:
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
        multiome_mudata_fname: pathlib.Path,
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
    mdata = mudata.read(multiome_mudata_fname.__str__())
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

class DEM(MotifEnrichmentResult):
    def __init__(
        self,
        foreground_regions: pr.PyRanges,
        background_regions: pr.PyRanges,
        name: str,
        species: Literal[
                "homo_sapiens", "mus_musculus", "drosophila_melanogaster"],
        adjpval_thr: float = 0.05,
        log2fc_thr: float = 1.0,
        mean_fg_thr: float = 0.0,
        motif_hit_thr: Optional[float] = None,
        path_to_motif_annotations: Optional[str] = None,
        annotation_version: str = 'v10nr_clust',
        annotation_to_use: list = ['Direct_annot', 'Motif_similarity_annot', 'Orthology_annot', 'Motif_similarity_and_Orthology_annot'],
        motif_similarity_fdr: float = 0.001,
        orthologous_identity_threshold: float = 0.0,
        motifs_to_use: Optional[list] = None):
        self.foreground_regions = foreground_regions
        self.background_regions = background_regions
        self.adjpval_thr = adjpval_thr
        self.log2fc_thr = log2fc_thr
        self.mean_fg_thr = mean_fg_thr
        self.motif_hit_thr = motif_hit_thr
        super().__init__(
            name = name,
            species = species,
            path_to_motif_annotations = path_to_motif_annotations,
            annotation_version = annotation_version,
            annotation_to_use = annotation_to_use,
            motif_similarity_fdr = motif_similarity_fdr,
            orthologous_identity_threshold = orthologous_identity_threshold,
            motifs_to_use = motifs_to_use)
    
    def run(self, dem_db: DEMDatabase):
        # Create logger
        level    = logging.INFO
        format   = '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'
        handlers = [logging.StreamHandler(stream=sys.stdout)]
        logging.basicConfig(level = level, format = format, handlers = handlers)
        log = logging.getLogger('DEM')

        log.info(f"Running DEM for {self.name}")
        foreground_region_to_db, foreground_scores = dem_db.get_scores(
            self.foreground_regions)
        background_region_to_db, background_scores = dem_db.get_scores(
            self.background_regions)
        self.regions_to_db = pd.concat(
            [foreground_region_to_db, background_region_to_db]) \
            .drop_duplicates() \
            .reset_index(drop = True)
        motif_names = foreground_scores.index
        background_scores = background_scores.loc[motif_names]
        foreground_scores_arr = np.array(foreground_scores)
        background_scores_arr = np.array(background_scores)
        
        # Perform Wilcoxon rank sum test
        stats, pvalues = ranksums_numba_multiple(
            X = foreground_scores_arr,
            Y = background_scores_arr)
        
        # Calculate log2FC
        logFC = get_log2_fc(
            fg_mat = foreground_scores_arr,
            bg_mat = background_scores_arr
        )
        # pvalue correction
        pvalues_adj = p_adjust_bh(pvalues)

        # Create result dataframe
        result = pd.DataFrame(
            data = {
                "Log2FC": logFC,
                "Adjusted_pval": pvalues_adj,
                "Mean_fg": mean_axis1(foreground_scores_arr),
                "Mean_bg": mean_axis1(background_scores_arr)},
            index = motif_names)

        if (result["Adjusted_pval"] <= self.adjpval_thr).sum() < 1:
            print("""
            WARNING: No significant pvalue, 
            we'll take the 1% top values to make the analysis go through
            """)
            self.adjpval_thr = min(
                self.adjpval_thr,
                result["Adjusted_pval"].quantile(0.01))

        # Threshold dataframe
        result = result.loc[
            np.logical_and.reduce(
                (
                    result["Adjusted_pval"] <= self.adjpval_thr,
                    result["Log2FC"]        >= self.log2fc_thr,
                    result["Mean_fg"]       >= self.mean_fg_thr
                )
            )
        ]

        result = result.sort_values([
            "Log2FC", "Adjusted_pval"], 
            ascending = [False, True])
        
        self.motif_enrichment = result
        log.info("Adding motif-to-TF annotation")
        self.add_motif_annotation()

        # Get motif hits
        result["Motif_hit_thr"] = None
        result["Motif_hits"] = None
        significant_motifs = result.index
        motif_hits = {}
        for motif in significant_motifs:
            if self.motif_hit_thr is None:
                thr = get_optimal_threshold_roc(
                    foreground_scores.loc[motif].values,
                    background_scores.loc[motif].values,)
            else:
                thr = self.motif_hit_thr 
            motif_hits[motif] = foreground_scores.loc[motif].loc[
                    foreground_scores.loc[motif] > thr].index.tolist()
            result.loc[motif, "Motif_hit_thr"] = thr
            result.loc[motif, "Motif_hits"] = len(motif_hits[motif])
        
        rs_motif_hits = {
            motif: list(set(
                self.regions_to_db.loc[
                    self.regions_to_db["Query"].isin(motif_hits[motif]),
                    "Target"].tolist()))
            for motif in motif_hits.keys()}
        self.motif_hits = {'Database': motif_hits, 'Region_set': rs_motif_hits}
        self.get_cistromes()

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
    dem_db_fname=dem_scores_fname,
    output_fname_dem_result=dem_results_path,
    output_fname_dem_html="",
    n_cpu=32,
    path_to_genome_annotation="aertslab/genomes/hg38/hg38_ensdb_v86.csv",
    temp_dir=temp_dir,
    species=species,
    fraction_overlap_w_dem_database=0.4,
    path_to_motif_annotations=path_to_motif_annotations,
    annotation_version=annotation_version,
    motif_similarity_fdr=0.05,
    orthologous_identity_threshold=0.0,
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
    cistarget_rankings_fname,
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
    temp_dir=temp_dir,
    species=species,
    annotations_to_use=["Direct_annot", "Orthology_annot"]
)

# Reformat cistromes results, merging DEM and Cistarget pairs
# Giving Direct and Extended cistromes
prepare_motif_enrichment_results(
    paths_to_motif_enrichment_results=[
        dem_results_path, cistarget_results_path],
    multiome_mudata_fname=mudata_path,
    out_file_direct_annotation=output_cistromes_annotations_direct,
    out_file_extended_annotation=output_cistromes_annotations_extended,
    out_file_tf_names=output_tf_names,
    direct_annotation=["Direct_annot"],
    extended_annotation=["Orthology_annot"])


import anndata as ad
import pandas as pd
import polars as pl

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
## Assert all regions are present only once
assert len(all_h5ad.var_names.unique())==len(all_h5ad.var_names)
all_h5ad.var = all_var

# Get max rankof motifs -- will give scenicplus TF-region score
from scenicplus.triplet_score import get_max_rank_of_motif_for_each_TF
max_rank = get_max_rank_of_motif_for_each_TF(all_h5ad, cistarget_rankings_fname)

# Keep only actual TF-region pairs
tfb_matrix = max_rank * all_h5ad.X.astype(bool)

# Stack and rename columns
tfb = tfb_matrix.stack().reset_index()
tfb.columns = ["cre", "tf", "score"]
tfb = tfb[tfb["score"]>0]

# scale rank in a significant pvalue -like range
min, max = tfb["score"].min(), tfb["score"].max()
tfb["score"] = (tfb["score"]-min)/(max - min)*(0.04-0.001)+0.001
tfb["score"] = -np.log10(tfb["score"])
tfb["cre"] = tfb["cre"].str.replace(":", "-")
# Save
tfb.to_csv(output_tfb)
