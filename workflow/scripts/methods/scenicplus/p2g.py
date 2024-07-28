import pathlib
from matplotlib import cm
import logging
import os
import subprocess
import sys
from typing import List, Literal, Optional, Tuple, Set, Union
import joblib
from matplotlib.colors import Normalize
from scipy.stats import pearsonr, spearmanr
from sklearn.ensemble import (ExtraTreesRegressor, GradientBoostingRegressor,
                              RandomForestRegressor)
from tqdm import tqdm

import joblib
import argparse
import scanpy as sc
import muon as mu
import pandas as pd
import polars as pl
import numpy as np
import scenicplus.data_wrangling.gene_search_space
# import scenicplus.enhancer_to_gene
# from scenicplus.enhancer_to_gene import calculate_regions_to_genes_relationships
import pyranges as pr


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

# Interact auto sql definition
INTERACT_AS = """table interact
"Interaction between two regions"
    (
    string chrom;      "Chromosome (or contig, scaffold, etc.). For interchromosomal, use 2 records"
    uint chromStart;   "Start position of lower region. For interchromosomal, set to chromStart of this region"
    uint chromEnd;     "End position of upper region. For interchromosomal, set to chromEnd of this region"
    string name;       "Name of item, for display.  Usually 'sourceName/targetName' or empty"
    uint score;        "Score from 0-1000."
    double value;      "Strength of interaction or other data value. Typically basis for score"
    string exp;        "Experiment name (metadata for filtering). Use . if not applicable"
    string color;      "Item color.  Specified as r,g,b or hexadecimal #RRGGBB or html color name, as in //www.w3.org/TR/css3-color/#html4."
    string sourceChrom;  "Chromosome of source region (directional) or lower region. For non-directional interchromosomal, chrom of this region."
    uint sourceStart;  "Start position source/lower/this region"
    uint sourceEnd;    "End position in chromosome of source/lower/this region"
    string sourceName;  "Identifier of source/lower/this region"
    string sourceStrand; "Orientation of source/lower/this region: + or -.  Use . if not applicable"
    string targetChrom; "Chromosome of target region (directional) or upper region. For non-directional interchromosomal, chrom of other region"
    uint targetStart;  "Start position in chromosome of target/upper/this region"
    uint targetEnd;    "End position in chromosome of target/upper/this region"
    string targetName; "Identifier of target/upper/this region"
    string targetStrand; "Orientation of target/upper/this region: + or -.  Use . if not applicable"
    )
"""

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


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--mudata', required=True)
parser.add_argument('-o', '--output', required=True)

parser.add_argument('-t', '--temp_dir', required=True)
parser.add_argument('-c', '--njobs', required=True, type=int)
parser.add_argument('-u', '--upstream', required=True)
parser.add_argument('-d', '--downstream', required=True)
parser.add_argument('-e', '--extend_tss', required=True)
parser.add_argument('-r', '--remove_promoters', required=True)
parser.add_argument('-z', '--use_gene_boundaries', required=True)

# regions_to_genes realtionships specific args
parser.add_argument('-p', '--importance_scoring_method', required=True)
parser.add_argument('-s', '--correlation_scoring_method', required=True)
parser.add_argument('-g', '--organism', required=True)
parser.add_argument('-m', '--chrom_sizes_m', required=True)
parser.add_argument('-j', '--chrom_sizes_h', required=True)

# upstream : Tuple[int, int]
#     Minimum and maximum (up until another gene) number of bps upstream of
#     TSS to include in the search space.
# downstream : Tuple[int, int]
#     Minimum and maximum (up until another gene) number of bps downstream of
#     gene's end to include in the search space.
# extend_tss : Tuple[int, int]
#     Number of bps to extend the TSS to define promoter region.
args = parser.parse_args()

organism = args.organism
if organism == 'hg38':
    chromsizes_fname = args.chrom_sizes_h
    species="hsapiens"
elif organism == 'mm10':
    chromsizes_fname = args.chrom_sizes_m
    species="mmusculus"
else:
    raise ValueError("Organism not hg38 nor mm10")

use_gene_boundaries = args.use_gene_boundaries
upstream = tuple([int(num) for num in args.upstream.split(' ')])
downstream = tuple([int(num) for num in args.downstream.split(' ')])
extend_tss = tuple([int(num) for num in args.extend_tss.split(' ')])
remove_promoters = args.remove_promoters
importance_scoring_method = args.importance_scoring_method
correlation_scoring_method = args.correlation_scoring_method
njobs = args.njobs

# mask_expr_dropout : bool
#     Whether to mask expression dropout.
# n_cpu : int
#     Number of parallel processes to run.
mask_expr_dropout = True

# Load data
mdata = mu.read_h5mu(args.mudata)

#Download chromosome sizes and gene body coordinates
result = scenicplus.data_wrangling.gene_search_space.download_gene_annotation_and_chromsizes(
        species=species,
        biomart_host="http://www.ensembl.org",
        use_ucsc_chromosome_style=True)

if type(result) is tuple:
    annot, chromsizes = result
else:
    annot = result
    print(
        "Chrosomome sizes was not found, please provide this information manually.")

# Calculate search space
mdata["atac"].var_names = mdata["atac"].var_names.str.split('-', 1).str[0] + ':' + mdata["atac"].var_names.str.split('-', 1).str[1]
mdata["rna"][:,mdata["rna"].X.sum(0)!=0]
mdata["atac"][:,mdata["atac"].X.sum(0)!=0]

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
        temp_dir=args.temp_dir,
        mask_expr_dropout=False,
        importance_scoring_method=importance_scoring_method,
        correlation_scoring_method=correlation_scoring_method,
        n_cpu=njobs,)

p2g = p2g[p2g["importance"]!=0]
p2g = p2g.loc[:, ["target", "region", "importance_x_rho"]]
p2g.columns = ["gene", "cre", "score"]
p2g = p2g[p2g["score"]!=0]
print(p2g)

p2g.to_csv(args.output)
