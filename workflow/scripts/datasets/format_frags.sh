#!/bin/bash

FILE_PATH=$1
SAMPLE_NAME=$(basename "$FILE_PATH" .frags.tsv.gz)
echo 'File path: ' $FILE_PATH
echo 'Sample id: ' $SAMPLE_NAME

# Process, modify, compress to bgzip format, and index
zcat "$FILE_PATH" | \
awk -v sample="$SAMPLE_NAME" '$0 !~ /^#/ {print $1"\t"$2"\t"$3"\t"sample"_"($4 ~ /-1$/ ? substr($4, 1, length($4)-2) : $4)"\t"$5}' | \
bgzip > "${FILE_PATH}_modified.frags.tsv.bgz"

# Index the bgzipped file with tabix
tabix -p bed "${FILE_PATH}_modified.frags.tsv.bgz"

# (Optional) Replace original file with the new bgzipped file
mv "${FILE_PATH}_modified.frags.tsv.bgz" "$FILE_PATH"
mv "${FILE_PATH}_modified.frags.tsv.bgz.tbi" "$FILE_PATH.tbi"
