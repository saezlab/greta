#!/bin/bash


# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --new_dir) new_dir="$2"; shift ;;
        --path_pre) path_pre="$2"; shift ;;
        --path_ann) path_ann="$2"; shift ;;
        --path_csz) path_csz="$2"; shift ;;
        --ext) ext="$2"; shift ;;
        --threads) threads="$2"; shift ;;
        --path_out) path_out="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

python -c "import mudata, sys; \
m = mudata.read(sys.argv[1]); \
m.mod['scRNA'] = m.mod['rna']; \
del m.mod['rna']; \
m.mod['scATAC'] = m.mod['atac']; \
del m.mod['atac']; \
m.mod['scATAC'].var_names = m.mod['scATAC'].var_names.str.replace('-', ':', 1); \
m.write(sys.argv[2])" $path_pre $new_dir/mdata.h5mu

scenicplus prepare_data search_spance \
--multiome_mudata_fname $new_dir/mdata.h5mu \
--gene_annotation_fname $path_ann \
--chromsizes_fname $path_csz \
--upstream 1000 $ext \
--downstream 1000 $ext \
--out_fname $new_dir/space.tsv

scenicplus grn_inference region_to_gene \
--multiome_mudata_fname $new_dir/mdata.h5mu \
--search_space_fname $new_dir/space.tsv \
--temp_dir $TMPDIR \
--out_region_to_gene_adjacencies $new_dir/rg_adj.tsv \
--n_cpu $threads

python -c "import pandas as pd; \
import sys; \
tab = pd.read_table(sys.argv[1]); \
tab = tab[tab['importance_x_rho'].abs() > 1e-16]; \
tab = tab[['region', 'target', 'importance_x_rho']]; \
tab['region'] = tab['region'].str.replace(':', '-'); \
tab.columns = ['cre', 'gene', 'score']; \
tab.to_csv(sys.argv[2], index=False)" $new_dir/rg_adj.tsv $path_out
