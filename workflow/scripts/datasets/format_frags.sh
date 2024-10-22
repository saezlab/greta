#!/bin/bash

FILE_PATH=$1
SAMPLE_NAME=$(basename "$FILE_PATH" .frags.tsv.gz)
zcat "$FILE_PATH" | \
awk -v sample="$SAMPLE_NAME" '$0 !~ /^#/ {print $1"\t"$2"\t"$3"\t"sample"_"($4 ~ /-1$/ ? substr($4, 1, length($4)-2) : $4)"\t"$5}' | \
gzip > "${FILE_PATH%.gz}_modified.frags.tsv.gz"
mv ${FILE_PATH%.gz}_modified.frags.tsv.gz $FILE_PATH
