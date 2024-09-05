import kiwisolver
import matplotlib.pyplot as plt
from matplotlib import cm
from typing import List, Literal, Optional, Tuple, Set, Dict, Union
import argparse
import pathlib
import numpy as np
import os
os.environ["NUMBA_CACHE_DIR"] = "/tmp"
import joblib
import muon as mu
import mudata
import pandas as pd
import logging
from scenicplus.cli.commands import _format_egrns

log = logging.getLogger("SCENIC+")

parser = argparse.ArgumentParser()
# input
parser.add_argument('-i', '--multiome_mudata_path', required=True)
parser.add_argument('-t', '--tfb_path', required=True)
parser.add_argument('-p', '--p2g_path', required=True)
parser.add_argument('-l', '--cistarget_db_path_human', required=True)
parser.add_argument('-k', '--cistarget_db_path_mouse', required=True)
# params
parser.add_argument('-m', '--method', required=True)
parser.add_argument('-c', '--n_cpu', required=True, type=int)
parser.add_argument('-d', '--temp_dir', required=True)
parser.add_argument('-g', '--organism', required=True)
parser.add_argument('--order_regions_to_genes_by', type=str, required=False,
        default="importance")
parser.add_argument('--order_TFs_to_genes_by', type=str, required=False,
        default="importance")
parser.add_argument('--gsea_n_perm', type=int, required=False,
        default=1000)
parser.add_argument('--quantile_thresholds_region_to_gene', type=float, required=False,
        nargs="*", default=[0.85, 0.90, 0.95])
parser.add_argument('--top_n_regionTogenes_per_gene', type=int, required=False,
        nargs="*", default=[5, 10, 15])
parser.add_argument('--top_n_regionTogenes_per_region', type=int, required=False,
        nargs="*", default=[])
parser.add_argument('--min_regions_per_gene', required=True, type=int)
parser.add_argument('--rho_threshold', type=float, required=False,
        default=0.05,)
parser.add_argument('--min_target_genes', type=int, required=False,
        default=10,)
# output
parser.add_argument('-s', '--tf_to_gene_prior_path', required=True)
parser.add_argument('-o', '--mdl_path', required=True)
parser.add_argument('-j', '--eRegulon_out_fname', required=True)

args = parser.parse_args()

# Parameters
multiome_mudata_path = pathlib.Path(args.multiome_mudata_path)
tfb_path = pathlib.Path(args.tfb_path)
p2g_path = pathlib.Path(args.p2g_path)
#tf_names_path = pathlib.Path(args.tf_names_path)
method = args.method
n_cpu = args.n_cpu
tf_to_gene_prior_path = pathlib.Path(args.tf_to_gene_prior_path)
temp_dir = pathlib.Path(args.temp_dir)
eRegulon_out_fname = pathlib.Path(args.eRegulon_out_fname)
mdl_path = pathlib.Path(args.mdl_path)
seed = 1

if args.organism == "hg38":
    ranking_db_fname = pathlib.Path(args.cistarget_db_path_human)
elif args.organism == "mm10":
    ranking_db_fname = pathlib.Path(args.cistarget_db_path_mouse)
else:
    raise ValueError("Invalid organism")

def infer_TF_to_gene(
        multiome_mudata_fname: pathlib.Path,
        tf_names: pathlib.Path,
        temp_dir: pathlib.Path,
        adj_out_fname: pathlib.Path,
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
    log.info(f"Saving TF to gene adjacencies to: {adj_out_fname.__str__()}")

    return adj
#    adj.to_csv(
 #       adj_out_fname,
  #      sep="\t", header = True, index = False)


def infer_grn(
        TF_to_gene_adj: pathlib.Path,
        region_to_gene_adj: pathlib.Path,
        cistromes: pathlib.Path,
        eRegulon_out_fname: pathlib.Path,
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
    from scenicplus.triplet_score import calculate_triplet_score
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
#    eRegulon_metadata = calculate_triplet_score(
 #       cistromes=cistromes,
  #      eRegulon_metadata=eRegulon_metadata,
   #     ranking_db_fname=ranking_db_fname)

    log.info(f"Saving network to {eRegulon_out_fname.__str__()}")
    return eRegulon_metadata #.to_csv(
        #eRegulon_out_fname,
        #sep="\t", header=True, index=False)

from scenicplus.utils import p_adjust_bh
from scenicplus.grn_builder.modules import (
    create_emodules, eRegulon, merge_emodules, RHO_THRESHOLD, TARGET_GENE_NAME)
from scenicplus.grn_builder.gsea import run_gsea
from scenicplus.grn_builder.gsea_approach import _run_gsea_for_e_module
import anndata
from tqdm import tqdm

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


# Load data
multiome_mudata = mu.read(multiome_mudata_path)
tfb = pd.read_csv(tfb_path)
p2g = pd.read_csv(p2g_path)

tf_names = list(tfb["tf"].unique())

# TF to gene relationships
tf_to_gene_prior = infer_TF_to_gene(
        multiome_mudata_fname=multiome_mudata_path,
        tf_names=tf_names,
        temp_dir=temp_dir,
        adj_out_fname=tf_to_gene_prior_path,
        method=method,
        n_cpu=n_cpu,
        seed=seed)

tf_to_gene_prior.to_csv("a.csv")

# Format p2g
p2g = p2g.rename(columns={
    "cre": "region",
    "gene": "target",
    "score": "importance_x_rho",
    })
p2g["importance_x_rho"][:100] = p2g["importance_x_rho"][:100]*(-1)
p2g["rho"] = np.sign(p2g["importance_x_rho"])
p2g["importance_x_abs_rho"] = np.abs(p2g["importance_x_rho"])
p2g["importance"] = p2g["importance_x_rho"]

p2g = p2g[["region", "target", "importance", "rho", "importance_x_rho", "importance_x_abs_rho"]]


# Format tfb and build cistromes
import polars as pl
import anndata as ad
tfb = tfb[["cre", "tf", "score"]]
tfb["cre"] = tfb["cre"].str.replace(":", "-")

tfb.columns = ["region", "tf", "score"]
cistromes = tfb[["region", "tf", "score"]]
cistromes = pl.DataFrame(cistromes).pivot(index="region", columns="tf", values="score").to_pandas()
cistromes = cistromes.set_index("region").astype(bool)
cistromes = ad.AnnData(cistromes)
print(cistromes.var_names, cistromes.obs_names)

print(p2g.head())
print(tf_to_gene_prior.head())

# Infer eGRN
mdl = infer_grn(
        TF_to_gene_adj=tf_to_gene_prior,
        region_to_gene_adj=p2g,
        cistromes=cistromes,
        eRegulon_out_fname=eRegulon_out_fname,
        ranking_db_fname=ranking_db_fname,
        is_extended=True,
        temp_dir=temp_dir,
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

mdl
