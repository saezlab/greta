#!/bin/bash


# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --input_pre) input_pre="$2"; shift ;;
        --input_p2g) input_p2g="$2"; shift ;;
        --input_frags)
            shift
            while [[ "$#" -gt 0 && "$1" != --* ]]; do
                input_frags+=("$1")
                shift
            done
            continue
            ;;
        --input_motif) input_motif="$2"; shift ;;
        --input_genome) input_genome="$2"; shift ;;
        --output_d) output_d="$2"; shift ;;
        --output_out) output_out="$2"; shift ;;
        --threads) threads="$2"; shift ;;
        --use_p2g) use_p2g="$2"; shift ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done && \
if [ $(wc -l < $input_p2g) -eq 1 ]; then
    echo "cre,tf,score" > "$output_out"
    mkdir -p "$output_d"
    exit 0
fi && \
mkdir -p "$output_d" && \
python workflow/scripts/mth/dictys/extract_data.py \
--pre_path "$input_pre" \
--p2g_path "$input_p2g" \
--exp_path "$output_d/expr.tsv.gz" \
--pks_path "$output_d/peaks.bed" \
--use_p2g "$use_p2g" && \
for b_file in $output_d/barcodes_*; do
    ctype=$(basename "$b_file" | sed 's/barcodes_//; s/.txt//')
    bam_name="$output_d/reads_$ctype.bam"
    bai_name="$output_d/reads_$ctype.bai"
    foot_name="$output_d/foot_$ctype.tsv.gz"
    motif_name="$output_d/motif_$ctype.tsv.gz"
    well_name="$output_d/well_$ctype.tsv.gz"
    homer_name="$output_d/homer_$ctype.tsv.gz"
    bind_name="$output_d/bind_$ctype.tsv.gz"
    tfb_bed="$output_d/tfb_$ctype.bed"
    echo "Processing $ctype"
    python workflow/scripts/mth/dictys/frag_to_bam.py --fnames "${input_frags[@]}" --barcodes $b_file | \
    samtools view -b | samtools sort -o "$bam_name" && samtools index "$bam_name" "$bai_name"
    python3 -m dictys chromatin wellington "$bam_name" "$bai_name" "$output_d/peaks.bed" "$foot_name" --nth "$threads" && \
    python3 -m dictys chromatin homer "$foot_name" "$input_motif" "$input_genome" "$output_d/expr.tsv.gz" "$motif_name" "$well_name" "$homer_name" && \
    python3 -m dictys chromatin binding "$well_name" "$homer_name" "$bind_name" && \
    zcat "$bind_name" | awk 'BEGIN {FS="\t"; OFS="\t"} NR>1 {split($2, coords, ":"); split(coords[3], pos, "-"); print coords[1], coords[2], coords[3], $1, $3}' | \
    bedtools intersect -a "$output_d/peaks.bed" -b stdin -wa -wb | awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2, $3, $7, $8}' > "$tfb_bed"
done && \
python -c "import pandas as pd, glob, sys, os; \
df = [pd.read_csv(f, sep='\t', header=None) for f in glob.glob(os.path.join(sys.argv[1], 'tfb_*.bed'))]; \
df = pd.concat(df); \
df.columns = ['chr', 'str', 'end', 'tf', 'score']; \
df['cre'] = df['chr'] + '-' + df['str'].astype(str) + '-' + df['end'].astype(str); \
df = df[['cre', 'tf', 'score']]; \
df = df.groupby(['cre', 'tf'], as_index=False)['score'].mean(); \
df.to_csv(sys.argv[2], index=False)" "$output_d" "$output_out"
