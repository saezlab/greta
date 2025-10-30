
rule get_mbrain_frags:
    threads: 2
    singularity: 'workflow/envs/gretabench.sif'
    input:
        img='workflow/envs/gretabench.sif',
    output:
        frag=temp(local(expand('dts/mm10/mbrain/{sample}.frags.tsv.gz', sample=config['dts']['mbrain']['samples']))),
    params:
        url_frags   = config['dts']['mbrain']['url']['atac_frags'],
        samples     = config['dts']['mbrain']['samples'],
    shell:
        """
        path_late="dts/mm10/mbrain"
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

rule download_mbrain:
    threads: 2
    singularity: 'workflow/envs/gretabench.sif'
    input:
        img='workflow/envs/gretabench.sif',
        gid=rules.gen_gid_ensmbl.output,
        frag=rules.get_mbrain_frags.output.frag,
    output:
        ann=temp(local('dts/mm10/mbrain/annot.csv')),
        gene_map=temp(local('dts/mm10/mbrain/mm10_gene_map.tsv')),
        rna=temp(local('dts/mm10/mbrain/rna.h5ad')),
    params:
        url_peaks   = config['dts']['mbrain']['url']['atac_peaks'],
        url_frags   = config['dts']['mbrain']['url']['atac_frags'],
        url_counts  = config['dts']['mbrain']['url']['atac_counts'],
        url_barc    = config['dts']['mbrain']['url']['barcodes'],
        url_ct      = config['dts']['mbrain']['url']['celltypes'],
        url_rna     = config['dts']['mbrain']['url']['rna_counts'],
        samples     = config['dts']['mbrain']['samples'],
    shell:
        """
        set -euo pipefail

        path_late="dts/mm10/mbrain"
        mkdir -p "$path_late"

        echo "[INFO] Downloading barcodes, celltypes, and RNA counts..."
        curl -L --fail --retry 3 --retry-delay 3 "{params.url_barc}" -o "$path_late/mbrain.barcodes.txt.gz"
        curl -L --fail --retry 3 --retry-delay 3 "{params.url_ct}"   -o "$path_late/mbrain.celltype.txt.gz"
        curl -L --fail --retry 3 --retry-delay 3 "{params.url_rna}"  -o "$path_late/rna.counts.txt.gz"

        echo "[INFO] Downloading ATAC peaks and counts..."
        curl -L --fail --retry 3 --retry-delay 3 "{params.url_peaks}" -o "$path_late/mbrain.peaks.bed.gz"
        curl -L --fail --retry 3 --retry-delay 3 "{params.url_counts}" -o "$path_late/mbrain.atac.counts.txt.gz"

        
        echo "[INFO] Preparing annot.csv..."
        gunzip -c "$path_late/mbrain.barcodes.txt.gz" > "$path_late/mbrain.barcodes.txt"
        gunzip -c "$path_late/mbrain.celltype.txt.gz" > "$path_late/mbrain.celltype.txt"

        sample_name=$(echo {params.samples} | tr -d "[],'" | xargs | cut -d' ' -f1)

        python workflow/scripts/dts/mbrain/annot.py \
          dts/mm10/mbrain/mbrain.barcodes.txt \
          dts/mm10/mbrain/mbrain.celltype.txt \
          "$sample_name"

        echo "[INFO] Building mm10 gene map (BiomaRt export)..."
        Rscript workflow/scripts/dts/mbrain/mouse_gene_map.R \
          dts/mm10/mbrain/mm10_gene_map.tsv

        echo "[INFO] Converting RNA counts to h5ad..."
        python workflow/scripts/dts/mbrain/rna_to_h5ad.py \
          dts/mm10/mbrain/rna.counts.txt.gz \
          dts/mm10/mbrain/annot.csv \
          dts/mm10/mbrain/mm10_gene_map.tsv \
          dts/mm10/mbrain/mbrain.celltype.txt \
          dts/mm10/mbrain/rna.h5ad

        echo "[CLEANUP] Removing intermediate downloads no longer needed..."
        rm -f "$path_late/mbrain.barcodes.txt.gz" \
              "$path_late/mbrain.barcodes.txt" \
              "$path_late/mbrain.celltype.txt" \
              "$path_late/rna.counts.txt.gz" \
              "$path_late/unmapped_genes.txt"
        """
    


rule callpeaks_mbrain:
    threads: 8
    singularity: 'workflow/envs/gretabench.sif'
    input:
        frags=rules.get_mbrain_frags.output.frag,
        annot=rules.download_mbrain.output.ann,
    output: peaks=temp(local('dts/mm10/mbrain/peaks.h5ad'))
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


rule annotate_mbrain:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        rna=rules.download_mbrain.output.rna,
        peaks=rules.callpeaks_mbrain.output.peaks,
        annot=rules.download_mbrain.output.ann,
    output: out='dts/mm10/mbrain/annotated.h5mu'
    shell:
        """
        python workflow/scripts/dts/mbrain/mbrain.py \
        -a {input.rna} \
        -b {input.peaks} \
        -c {input.annot} \
        -d {output.out}
        """
