import argparse
import os
os.environ['NUMBA_CACHE_DIR'] = '/tmp/'
import scanpy as sc
import muon as mu
import pandas as pd
import scenicplus

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--mudata', required=True)
parser.add_argument('-o', '--output', required=True)
parser.add_argument('-a', '--annot', required=True)
parser.add_argument('-n', '--chromsizes', required=True)
parser.add_argument('-t', '--temp_dir', required=True)
parser.add_argument('-c', '--njobs', required=True)
parser.add_argument('-u', '--upstream', required=True)
parser.add_argument('-d', '--downstream', required=True)
parser.add_argument('-e', '--extend_tss', required=True)
parser.add_argument('-r', '--remove_promoters', required=True)
parser.add_argument('-g', '--use_gene_boundaries', required=True)

# regions_to_genes realtionships specific args
parser.add_argument('-p', '--importance_scoring_method', required=True)
parser.add_argument('-s', '--correlation_scoring_method', required=True)
parser.add_argument('-m', '--chrom_sizes_m', required=True)
parser.add_argument('-h', '--chrom_sizes_h', required=True)
parser.add_argument('-a', '--annot_human', required=True)
parser.add_argument('-b', '--annot_mouse', required=True)

# upstream : Tuple[int, int]
#     Minimum and maximum (up until another gene) number of bps upstream of
#     TSS to include in the search space.
# downstream : Tuple[int, int]
#     Minimum and maximum (up until another gene) number of bps downstream of
#     gene's end to include in the search space.
# extend_tss : Tuple[int, int]
#     Number of bps to extend the TSS to define promoter region.
args = parser.parse_args()

organism = args['organism']
if organism == 'hsapiens':
    chromsizes_fname = args['chrom_sizes_h']
    annot_fname = args['annot_human']
if organism == 'mmusculus':
    chromsizes_fname = args['chrom_sizes_m']
    annot_fname = args['annot_mouse']

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
gene_annotation = pd.read_table(annot_fname)
chromsizes = pd.read_table(chromsizes_fname)


# Calculate search space
search_space = scenicplus.data_wrangling.gene_search_space.get_search_space(
        scplus_region=set(mdata["scATAC"].var_names),
        scplus_genes=set(mdata["scRNA"].var_names),
        gene_annotation=gene_annotation,
        chromsizes=chromsizes,
        use_gene_boundaries=use_gene_boundaries,
        upstream=upstream,
        downstream=downstream,
        extend_tss=extend_tss,
        remove_promoters=remove_promoters)

# Calculate regions to genes relationships
scenicplus.enhancer_to_gene.calculate_regions_to_genes_relationships(
        df_exp_mtx=mdata["scRNA"].to_df(),
        df_acc_mtx=mdata["scATAC"].to_df(),
        search_space=search_space,
        temp_dir=args.temp_dir,
        mask_expr_dropout=mask_expr_dropout,
        importance_scoring_method=importance_scoring_method,
        correlation_scoring_method=correlation_scoring_method,
        n_cpu=njobs,)
