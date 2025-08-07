#!/usr/bin/env bash
set -euo pipefail

echo "Processing KPMP zips..."

for zip in dts/KPMP/KPMP_*.zip; do
    echo "Trying to unzip: $zip"
    python3 -c "import zipfile; zipfile.ZipFile('$zip').extractall('dts/KPMP/')"
    echo "Unzip worked: $zip"

    filename=$(basename "$zip")
    folder="dts/KPMP/$(echo "$filename" | sed 's/\.zip$//')"
    id=$(echo "$folder" | cut -d'_' -f2)

    echo "Folder: $folder"
    echo "ID: $id"
    echo "${folder}/atac_fragments.tsv.gz"

    gzip -cd "${folder}/atac_fragments.tsv.gz" | \
    awk -v ID="$id" 'BEGIN {OFS="\t"}
        $0 !~ /^#/ {
            gsub("-1$", "", $4)
            $4 = ID "_" $4
            print
        }' > "dts/KPMP/${id}.frags.tsv"
done

for f in dts/KPMP/*.frags.tsv; do
    bgzip -f -c "$f" > "$f.gz"
    tabix -p bed "$f.gz"
    rm "$f"
    rm "$f.gz.tbi"
done

rm dts/KPMP/KPMP_*.zip
rm -r dts/KPMP/KPMP_*/