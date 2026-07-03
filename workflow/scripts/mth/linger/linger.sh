#!/bin/bash

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --linger_GRN) linger_GRN="$2"; shift ;;
        --out_dir) out_dir="$2"; shift ;;
        --path_mdata) path_mdata="$2"; shift ;;
        --version) version="$2"; shift ;;
        --genome) genome="$2"; shift ;;
        --mode) mode="$2"; shift ;;
        --path_out) path_out="$2"; shift ;;
        *) echo echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

python workflow/scripts/mth/linger/linger.py \
-g $linger_GRN \
-d $out_dir \
-m $path_mdata \
-v $version \
-gen $genome \
-md $mode

python workflow/scripts/mth/linger/grn.py \
-d $out_dir \
-o $path_out