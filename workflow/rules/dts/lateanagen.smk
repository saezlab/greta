
rule get_lateanagen_frags:
    threads: 2
    singularity: 'workflow/envs/gretabench.sif'
    input:
        img='workflow/envs/gretabench.sif',
    output:
        frag=temp(local(expand('dts/lateanagen/{sample}.frags.tsv.gz', sample=config['dts']['lateanagen']['samples']))),
    params:
        url_frags   = config['dts']['lateanagen']['url']['atac_frags'],
        samples     = config['dts']['lateanagen']['samples'],
    shell:
        """
        path_late="dts/lateanagen"
        samples=$(echo {params.samples} | tr -d "[],'")
        for s in $samples; do
            echo "[INFO] Downloading fragment file for $s ..."
            curl -L --fail --retry 3 --retry-delay 3 "{params.url_frags}" \
                -o "$path_late/${{s}}.frags.bed.gz"
            
            echo "[INFO] Converting BED format fragments to TSV format and sorting..."
            # The fragment file is in BED format with columns: chr, start, end, barcode, count
            # Need to sort by barcode and position for snapatac2
            gunzip -c "$path_late/${{s}}.frags.bed.gz" | \
                sort -k4,4 -k1,1 -k2,2n | \
                gzip > "$path_late/${{s}}.frags.tsv.gz"
            
            rm "$path_late/${{s}}.frags.bed.gz"
        done
        """

rule download_lateanagen:
    threads: 2
    singularity: 'workflow/envs/gretabench.sif'
    input:
        img='workflow/envs/gretabench.sif',
        gid=rules.gen_gid_ensmbl.output,
        frag=rules.get_lateanagen_frags.output.frag,
    output:
        ann=temp(local('dts/lateanagen/annot.csv')),
        gene_map=temp(local('dts/lateanagen/mm10_gene_map.tsv')),
        rna=temp(local('dts/lateanagen/rna.h5ad')),
    params:
        url_peaks   = config['dts']['lateanagen']['url']['atac_peaks'],
        url_frags   = config['dts']['lateanagen']['url']['atac_frags'],
        url_counts  = config['dts']['lateanagen']['url']['atac_counts'],
        url_barc    = config['dts']['lateanagen']['url']['barcodes'],
        url_ct      = config['dts']['lateanagen']['url']['celltypes'],
        url_rna     = config['dts']['lateanagen']['url']['rna_counts'],
        samples     = config['dts']['lateanagen']['samples'],
    shell:
        """
        set -euo pipefail

        path_late="dts/lateanagen"
        mkdir -p "$path_late"

        echo "[INFO] Downloading barcodes, celltypes, and RNA counts..."
        curl -L --fail --retry 3 --retry-delay 3 "{params.url_barc}" -o "$path_late/lateanagen.barcodes.txt.gz"
        curl -L --fail --retry 3 --retry-delay 3 "{params.url_ct}"   -o "$path_late/lateanagen.celltype.txt.gz"
        curl -L --fail --retry 3 --retry-delay 3 "{params.url_rna}"  -o "$path_late/rna.counts.txt.gz"

        echo "[INFO] Downloading ATAC peaks and counts..."
        curl -L --fail --retry 3 --retry-delay 3 "{params.url_peaks}" -o "$path_late/lateanagen.peaks.bed.gz"
        curl -L --fail --retry 3 --retry-delay 3 "{params.url_counts}" -o "$path_late/lateanagen.atac.counts.txt.gz"

        
        echo "[INFO] Preparing annot.csv..."
        gunzip -c "$path_late/lateanagen.barcodes.txt.gz" > "$path_late/lateanagen.barcodes.txt"
        gunzip -c "$path_late/lateanagen.celltype.txt.gz" > "$path_late/lateanagen.celltype.txt"

        sample_name=$(echo {params.samples} | tr -d "[],'" | xargs | cut -d' ' -f1)

        python workflow/scripts/dts/lateanagen/annot.py \
          dts/lateanagen/lateanagen.barcodes.txt \
          dts/lateanagen/lateanagen.celltype.txt \
          "$sample_name"

        echo "[INFO] Building mm10 gene map (BiomaRt export)..."
        Rscript workflow/scripts/dts/lateanagen/mouse_gene_map.R \
          dts/lateanagen/mm10_gene_map.tsv

        echo "[INFO] Converting RNA counts to h5ad..."
        python workflow/scripts/dts/lateanagen/rna_to_h5ad.py \
          dts/lateanagen/rna.counts.txt.gz \
          dts/lateanagen/annot.csv \
          dts/lateanagen/mm10_gene_map.tsv \
          dts/lateanagen/lateanagen.celltype.txt \
          dts/lateanagen/rna.h5ad

        echo "[CLEANUP] Removing intermediate downloads no longer needed..."
        rm -f "$path_late/lateanagen.barcodes.txt.gz" \
              "$path_late/lateanagen.barcodes.txt" \
              "$path_late/lateanagen.celltype.txt" \
              "$path_late/rna.counts.txt.gz" \
              "$path_late/unmapped_genes.txt"
        """
    


rule callpeaks_lateanagen:
    threads: 8
    singularity: 'workflow/envs/gretabench.sif'
    input:
        frags=rules.get_lateanagen_frags.output.frag,
        annot=rules.download_lateanagen.output.ann,
    output: peaks=temp(local('dts/lateanagen/peaks.h5ad'))
    resources:
        mem_mb=64000,
        runtime=360,
    shell:
        """
        python workflow/scripts/dts/callpeaks_mm10.py \
        -f {input.frags} \
        -a {input.annot} \
        -t $TMPDIR \
        -n {threads} \
        -o {output.peaks}
        """


rule annotate_lateanagen:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        rna=rules.download_lateanagen.output.rna,
        peaks=rules.callpeaks_lateanagen.output.peaks,
        annot=rules.download_lateanagen.output.ann,
    output: out='dts/lateanagen/annotated.h5mu'
    shell:
        """
        python workflow/scripts/dts/lateanagen/lateanagen.py \
        -a {input.rna} \
        -b {input.peaks} \
        -c {input.annot} \
        -d {output.out}
        """
