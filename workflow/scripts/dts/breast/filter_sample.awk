BEGIN { OFS="\t" }
$1 !~ /^chr/ { next }
{
    split($4, arr, "_");          # arr[1] = barcode, arr[2] = pool id
    gsub("-1", "", arr[1]);       # remove trailing -1
    current_pool = "ID" arr[2];

    # Only process lines matching the pool passed from parallel
    if (current_pool == pool) {
        new_barcode = current_pool "_" arr[1];
        print $1, $2, $3, new_barcode, $5 >> output_file
        fflush(output_file)      # optional, ensures immediate write
    }
}