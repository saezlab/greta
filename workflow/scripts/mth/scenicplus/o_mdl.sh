#!/bin/bash


# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --new_dir) new_dir="$2"; shift ;;
        --path_mdata) path_mdata="$2"; shift ;;
        --path_blist) path_blist="$2"; shift ;;
        --ntopics) ntopics="$2"; shift ;;
        --path_ann) path_ann="$2"; shift ;;
        --path_csz) path_csz="$2"; shift ;;
        --ext) ext="$2"; shift ;;
        --path_rnk) path_rnk="$2"; shift ;;
        --path_man) path_man="$2"; shift ;;
        --path_scr) path_scr="$2"; shift ;;
        --threads) threads="$2"; shift ;;
        --path_out) path_out="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# pre
python workflow/scripts/mth/scenicplus/topics.py \
-m $path_mdata \
-b $path_blist \
-n $ntopics \
-c $threads \
-d $new_dir \
-o $new_dir/mdata.h5mu

# p2g
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

# tfb
pycistopic topic_modeling mallet binarize \
--target region \
--method otsu \
--smooth True \
--nbins 100 \
--regions $new_dir/cres.txt \
--n_topics $ntopics \
-o $new_dir/topics/topic_ \
-p $new_dir/topics/otsu

scenicplus grn_inference motif_enrichment_cistarget \
--region_set_folder $new_dir/topics/ \
--cistarget_db_fname $path_rnk \
--output_fname_cistarget_result $new_dir/cistarget.hdf5 \
--path_to_motif_annotations $path_man \
--annotations_to_use Direct_annot Orthology_annot \
--temp_dir $TMPDIR \
--species homo_sapiens \
--n_cpu 4

scenicplus grn_inference motif_enrichment_dem  \
--region_set_folder $new_dir/topics/ \
--dem_db_fname $path_scr \
--output_fname_dem_result $new_dir/dem.hdf5 \
--max_bg_regions 500 \
--balance_number_of_promoters \
--genome_annotation $path_ann \
--motif_hit_thr 3 \
--path_to_motif_annotations $path_man \
--annotations_to_use Direct_annot Orthology_annot \
--temp_dir $TMPDIR \
--species homo_sapiens \
--n_cpu 4

scenicplus prepare_data prepare_menr \
--paths_to_motif_enrichment_results $new_dir/cistarget.hdf5 $new_dir/dem.hdf5 \
--multiome_mudata_fname $new_dir/mdata.h5mu \
--out_file_tf_names $new_dir/tfs.txt \
--out_file_direct_annotation $new_dir/direct.h5ad \
--out_file_extended_annotation $new_dir/extended.h5ad \
--direct_annotation Direct_annot Orthology_annot \
--extended_annotation Orthology_annot

scenicplus grn_inference TF_to_gene \
--multiome_mudata_fname $new_dir/mdata.h5mu \
--tf_names $new_dir/tfs.txt \
--temp_dir $TMPDIR \
--out_tf_to_gene_adjacencies $new_dir/tg_adj.tsv \
--method GBM \
--n_cpu $threads

scenicplus grn_inference eGRN \
--TF_to_gene_adj_fname $new_dir/tg_adj.tsv \
--region_to_gene_adj_fname $new_dir/rg_adj.tsv \
--cistromes_fname $new_dir/direct.h5ad \
--ranking_db_fname $path_rnk \
--eRegulon_out_fname $new_dir/egrn.tsv \
--temp_dir $TMPDIR \
--min_target_genes 10 \
--n_cpu $threads

python workflow/scripts/mth/scenicplus/egrn.py $new_dir/egrn.tsv $path_out
