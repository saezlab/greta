#!/usr/bin/env bash
pool="$1"
path_breast="$2"
path_tmp="$3"
output_file="$path_breast/${pool}.frags.tsv"
gunzip -c "$path_tmp" | awk -v pool="$pool" -v output_file="$output_file" -f workflow/scripts/dts/breast/filter_sample.awk