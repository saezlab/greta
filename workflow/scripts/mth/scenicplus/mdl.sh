#!/bin/bash


# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --new_dir) new_dir="$2"; shift ;;
        --path_pre) path_pre="$2"; shift ;;
        --path_p2g) path_p2g="$2"; shift ;;
        --path_tfb) path_tfb="$2"; shift ;;
        --path_rnk) path_rnk="$2"; shift ;;
        --threads) threads="$2"; shift ;;
        --path_out) path_out="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

if [ $(wc -l < $path_p2g) -eq 1 ] || [ $(wc -l < $path_tfb) -eq 1 ]; then
    echo "source,target,score" > "$path_out"
    exit 0
fi

# Transform pre to scenicplus mdata format
python workflow/scripts/mth/scenicplus/mdata.py $path_pre $new_dir/mdata.h5mu

# Extract unique tfs
python -c "import pandas as pd; \
import sys; \
tfb = pd.DataFrame(pd.read_csv(sys.argv[1])['tf'].unique().reshape(-1, 1), columns=['tf']); \
tfb.to_csv(sys.argv[2], index=False, header=False)" $path_tfb $new_dir/tfs.txt

# Infer tg links
scenicplus grn_inference TF_to_gene \
--multiome_mudata_fname $new_dir/mdata.h5mu \
--tf_names $new_dir/tfs.txt \
--temp_dir $TMPDIR \
--out_tf_to_gene_adjacencies $new_dir/tg_adj.tsv \
--method GBM \
--n_cpu $threads

# Transform p2g to scenicplus format
python -c "import pandas as pd; \
import numpy as np; \
import sys; \
p2g = pd.read_csv(sys.argv[1]); \
p2g['region'] = p2g['cre'].str.replace('-', ':', 1); \
p2g['target'] = p2g['gene']; \
p2g['importance'] = p2g['score'].abs(); \
p2g['rho'] = np.sign(p2g['score']); \
p2g['importance_x_rho'] = p2g['score']; \
p2g['importance_x_abs_rho'] = p2g['score']; \
p2g = p2g[['region', 'target', 'importance', 'rho', 'importance_x_rho', 'importance_x_abs_rho']]; \
p2g.to_csv(sys.argv[2], sep='\\t', index=False)" $path_p2g $new_dir/rg_adj.tsv

# Transform tfb to scenicplus format
python workflow/scripts/mth/scenicplus/motifs.py $path_tfb $path_rnk $new_dir/motifs.h5ad

# Egrn inference
dichotomize=$( python -c "import sys, pandas; print('') if (pandas.read_csv(sys.argv[1])['score'] < 0).any() else print('--do_not_rho_dichotomize_r2g --do_not_rho_dichotomize_eRegulon');" $path_p2g )
echo "$dichotomize"
scenicplus grn_inference eGRN \
--TF_to_gene_adj_fname $new_dir/tg_adj.tsv \
--region_to_gene_adj_fname $new_dir/rg_adj.tsv \
--cistromes_fname $new_dir/motifs.h5ad \
--ranking_db_fname $path_rnk \
--eRegulon_out_fname $new_dir/egrn.tsv \
--temp_dir $TMPDIR \
--min_target_genes 10 \
--n_cpu $threads $dichotomize

# Transform grn into greta format
python workflow/scripts/mth/scenicplus/egrn.py $new_dir/egrn.tsv $path_out
