import yaml
import sys
import os
import warnings
#os.environ['NUMBA_CACHE_DIR'] = '/tmp/'
import pathlib
import argparse
import scanpy as sc
import anndata as ad
import pandas as pd
import muon as mu
import anndata
# Init args
parser = argparse.ArgumentParser()
    # Data and folders
parser.add_argument('-f', '--frags',  nargs='+', help='Path to fragments file')
parser.add_argument('-i', '--mudata', type=str, help='Path to metadata file')
parser.add_argument('-o', '--out', type=str, help='Path to output directory')
parser.add_argument('-t', '--temp_dir', type=str, help='Path to temp directory')
parser.add_argument('--ray_tmp_dir', type=str, help='Path to ray tmp directory')
parser.add_argument('--cistopic_path', type=str, help='Path to cisTopic object')
    # tfb intermediate files
parser.add_argument('--annotation_direct_path', type=str, help='Path to direct annotation')
parser.add_argument('--annotation_extended_path', type=str, help='Path to extended annotation')
parser.add_argument('--tf_names_path', type=str, help='Path to TF names')
parser.add_argument('--search_space_path', type=str, help='Search space path')
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
parser.add_argument('--output_config', type=str, required=True)
args = vars(parser.parse_args())

# Set up variables
frags = args['frags']
mudata_file = args['mudata']
output_fname = os.path.join('../../', args['out'])
temp_dir = os.path.join('../../', args['temp_dir'])
ray_tmp_dir = os.path.join("../..", args['ray_tmp_dir'])
njobs = args['njobs']
organism = args['organism']
dem_results_path = os.path.join('../../', args['dem_results_path'])
cistarget_results_path = os.path.join('../../', args['cistarget_results_path'])

if organism == 'hg38':
    organism = "hsapiens"
    species = "homo_sapiens",
    chromsizes_fname = os.path.join('../../', args['chrom_sizes_h'])
    annot_fname = os.path.join('../../', args['annot_human'])
    cistarget_ranking_db_fname = os.path.join('../../', args["cistarget_rankings_human"])
    cistarget_score_db_fname = os.path.join('../../', args["cistarget_scores_human"])
    path_to_motif_annotations = os.path.join('../../', args["path_to_motif_annotations_human"])
elif organism == 'mm10':
    organism = "mmusculus"
    chromsizes_fname = args['chrom_sizes_m']
    annot_fname = args['annot_mouse']
    cistarget_ranking_db_fname = str(pathlib.Path(args["cistarget_rankings_mouse"]))
    cistarget_score_db_fname = str(pathlib.Path(args["cistarget_scores_mouse"]))
    path_to_motif_annotations = args["path_to_motif_annotations_mouse"]
    species = "mus_musculus"
annotation_version = "v10nr_clust"


with open(args["output_config"]) as file:
    config_dic = yaml.load(file, Loader=yaml.BaseLoader)

print(config_dic)
cis_topic_tmp_dir = os.path.join(temp_dir, 'cisTopic')
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
GEX_anndata_fname = os.path.join(args["temp_dir"], 'GEX_anndata.h5ad')

# Save GEW anndata from mudata
mdata = mu.read_h5mu(mudata_file)
rna = mdata["rna"]

# Add a raw slot in the rna adata
raw_counts = rna.layers["counts"]  # Access the sparse raw counts
cell_names = rna.obs.index  # Cell names
gene_names = rna.var.index   # Gene names
raw_adata = anndata.AnnData(X=raw_counts, obs=rna.obs.copy(), var=rna.var.copy())
rna.raw = raw_adata

rna.write_h5ad(GEX_anndata_fname)

mudata_file= os.path.join("../..", mudata_file)

config_dic["input_data"]["cisTopic_obj_fname"] = os.path.join("../../", args["cistopic_path"])
config_dic["input_data"]["GEX_anndata_fname"] = os.path.join("../..", GEX_anndata_fname)
config_dic["input_data"]["region_set_folder"] = cis_topic_tmp_dir
config_dic["input_data"]["ctx_db_fname"] = cistarget_ranking_db_fname
config_dic["input_data"]["dem_db_fname"] = cistarget_score_db_fname
config_dic["input_data"]["path_to_motif_annotations"] = path_to_motif_annotations

config_dic["output_data"]["combined_GEX_ACC_mudata"] = os.path.join(temp_dir, 'combined_GEX_ACC_mudata.h5mu')
config_dic["output_data"]["dem_result_fname"] = dem_results_path
config_dic["output_data"]["ctx_result_fname"] = cistarget_results_path
config_dic["output_data"]["output_fname_dem_html"] = os.path.join(temp_dir, 'o_scenciplus_dem_results.html')
config_dic["output_data"]["output_fname_ctx_html"] = os.path.join(temp_dir, 'o_scenicplus_ctx_results.html')
config_dic["output_data"]["cistromes_direct"] = args['annotation_direct_path']
config_dic["output_data"]["cistromes_extended"] = args['annotation_extended_path']
config_dic["output_data"]["tf_names"] = args['tf_names_path']
config_dic["output_data"]["genome_annotation"] = annot_fname
config_dic["output_data"]["chromsizes"] = chromsizes_fname
config_dic["output_data"]["search_space"] = args['search_space_path']
config_dic["output_data"]["tf_to_gene_adjacencies"] = os.path.join(temp_dir, 'o_scenciplus_tf_to_gene_adjacencies.tsv')
config_dic["output_data"]["region_to_gene_adjacencies"] = os.path.join(temp_dir, 'o_scenciplus_region_to_gene_adjacencies.tsv')
config_dic["output_data"]["eRegulons_direct"] = os.path.join(temp_dir, 'o_scenciplus_eRegulons_direct.tsv')
config_dic["output_data"]["eRegulons_extended"] = os.path.join(temp_dir, 'o_scenciplus_eRegulons_extended.tsv')
config_dic["output_data"]["AUCell_direct"] = os.path.join(temp_dir, 'o_scenciplus_AUCell_direct.h5mu')
config_dic["output_data"]["AUCell_extended"] = os.path.join(temp_dir, 'o_scenciplus_AUCell_extended.h5mu')
config_dic["output_data"]["scplus_mdata"] = os.path.join(temp_dir, 'o_scenciplus_scplus_mdata.h5mu')

config_dic["params_general"]["n_cpu"] = njobs
config_dic["params_general"]["seed"] = 42
config_dic["params_general"]["temp_dir"] = args['temp_dir']

#config_dic["params_data_preparation"]["bc_transform_func"] = "lambda x: x"
config_dic["params_data_preparation"]["is_multiome"] = True
config_dic["params_data_preparation"]["key_to_group_by"] = ""
nr_cells_per_metacells = 10
config_dic["params_data_preparation"]["direct_annotation"] = "Direct_annot"
config_dic["params_data_preparation"]["extended_annotation"] = "Orthology_annot"
config_dic["params_data_preparation"]["species"] = organism
config_dic["params_data_preparation"]["biomart_host"] = "http://ensembl.org/"
config_dic["params_data_preparation"]["search_space_upstream"] = '{}'.format(args['search_space_upstream'])
config_dic["params_data_preparation"]["search_space_downstream"] =  '{}'.format(args['search_space_downstream'])
config_dic["params_data_preparation"]["search_space_extend_tss"] = '{}'.format(args['search_space_extend_tss'])

#config_dic["params_motif_enrichment"]["species"] = species
config_dic["params_motif_enrichment"]["annotation_version"] = annotation_version
config_dic["params_motif_enrichment"]["motif_similarity_fdr"] = 0.05
config_dic["params_motif_enrichment"]["orthologous_identity_threshold"] = 0.8
config_dic["params_motif_enrichment"]["annotations_to_use"] = "Direct_annot Orthology_annot"
config_dic["params_motif_enrichment"]["fraction_overlap_w_dem_database"] = 0.4
config_dic["params_motif_enrichment"]["dem_max_bg_regions"] = 500
config_dic["params_motif_enrichment"]["dem_balance_number_of_promoters"] = True
config_dic["params_motif_enrichment"]["dem_promoter_space"] = 1_000
config_dic["params_motif_enrichment"]["dem_adj_pval_thr"] = 0.05
config_dic["params_motif_enrichment"]["dem_log2fc_thr"] = 1.0
config_dic["params_motif_enrichment"]["dem_mean_fg_thr"] = 0.0
config_dic["params_motif_enrichment"]["dem_motif_hit_thr"] = 3.0
config_dic["params_motif_enrichment"]["fraction_overlap_w_ctx_database"] = 0.4
config_dic["params_motif_enrichment"]["ctx_auc_threshold"] = 0.005
config_dic["params_motif_enrichment"]["ctx_nes_threshold"] = 3.0
config_dic["params_motif_enrichment"]["ctx_rank_threshold"] = 0.05

config_dic["params_inference"]["tf_to_gene_importance_method"] = args['method_mdl']
config_dic["params_inference"]["region_to_gene_importance_method"] = args["region_to_gene_importance_method"]
config_dic["params_inference"]["region_to_gene_correlation_method"] = args["region_to_gene_correlation_method"]
config_dic["params_inference"]["order_regions_to_genes_by"] = args['order_regions_to_genes_by']
config_dic["params_inference"]["order_TFs_to_genes_by"] = args['order_TFs_to_genes_by']
config_dic["params_inference"]["gsea_n_perm"] = args['gsea_n_perm']
print(args["quantile_thresholds_region_to_gene"])
config_dic["params_inference"]["quantile_thresholds_region_to_gene"] = ' '.join([str(v) for v in args['quantile_thresholds_region_to_gene']])
config_dic["params_inference"]["top_n_regionTogenes_per_gene"] = ' '.join([str(v) for v in args['top_n_regionTogenes_per_gene']])
config_dic["params_inference"]["top_n_regionTogenes_per_region"] = ' '.join([str(v) for v in args['top_n_regionTogenes_per_region']])
config_dic["params_inference"]["min_regions_per_gene"] = args['min_regions_per_gene']
config_dic["params_inference"]["rho_threshold"] = args['rho_threshold']
config_dic["params_inference"]["min_target_genes"] = args['min_target_genes']


with open(args["output_config"], 'w') as file:
    yaml.dump(config_dic, file, default_flow_style=False)
