localrules: annotate_retina


rule extract_retina_rds:
    threads: 2
    singularity: 'workflow/envs/figr.sif'
    input:
        img='workflow/envs/figr.sif',
    output:
        annot=temp(local('dts/retina/annot.csv')),
        rna_mtx=temp(local('dts/retina/rna_counts.mtx')),
        rna_genes=temp(local('dts/retina/rna_genes.txt')),
        rna_barcodes=temp(local('dts/retina/rna_barcodes.txt')),
        peaks_mtx=temp(local('dts/retina/peaks_counts.mtx')),
        peak_names=temp(local('dts/retina/peak_names.txt')),
        peak_barcodes=temp(local('dts/retina/peak_barcodes.txt')),
    params:
        url_atac=config['dts']['retina']['url']['atac_rds'],
        url_rna=config['dts']['retina']['url']['rna_rds'],
        url_multiomics=config['dts']['retina']['url']['geneexpr_peaks'],
        sample=config['dts']['retina']['samples'][0],
    shell:
        """
        set -euo pipefail
        
        path_retina="dts/retina"
        mkdir -p "$path_retina"
        
        echo "[INFO] Downloading ATAC RDS file..."
        wget "{params.url_atac}" -O "$path_retina/atac.RDS.gz"
        gunzip "$path_retina/atac.RDS.gz"
        
        echo "[INFO] Downloading RNA RDS file..."
        wget "{params.url_rna}" -O "$path_retina/rna.RDS.gz"
        gunzip "$path_retina/rna.RDS.gz"
        
        echo "[INFO] Downloading multiomics RDS file (with peaks)..."
        wget "{params.url_multiomics}" -O "$path_retina/multiomics.RDS.gz"
        gunzip "$path_retina/multiomics.RDS.gz"
        
        echo "[INFO] Extracting data from RDS files..."
        Rscript workflow/scripts/dts/retina/extract_rds.R \
            "$path_retina/atac.RDS" \
            "$path_retina/rna.RDS" \
            "$path_retina/multiomics.RDS" \
            "$path_retina" \
            "{params.sample}"
        
        echo "[CLEANUP] Removing RDS files..."
        rm -f "$path_retina"/*.RDS
        
        echo "[INFO] RDS extraction complete!"
        """


rule download_retina:
    threads: 2
    singularity: 'workflow/envs/gretabench.sif'
    input:
        img='workflow/envs/gretabench.sif',
        gid='dbs/mm10/gen/gid/ensembl.csv',
        annot=rules.extract_retina_rds.output.annot,
        rna_mtx=rules.extract_retina_rds.output.rna_mtx,
        rna_genes=rules.extract_retina_rds.output.rna_genes,
        rna_barcodes=rules.extract_retina_rds.output.rna_barcodes,
        peaks_mtx=rules.extract_retina_rds.output.peaks_mtx,
        peak_names=rules.extract_retina_rds.output.peak_names,
        peak_barcodes=rules.extract_retina_rds.output.peak_barcodes,
    output:
        rna=temp(local('dts/retina/rna.h5ad')),
        peaks=temp(local('dts/retina/peaks.h5ad')),
    shell:
        """
        set -euo pipefail
        
        path_retina="dts/retina"
        
        echo "[INFO] Processing RNA counts to h5ad..."
        python workflow/scripts/dts/retina/process_rna.py \
            {input.rna_mtx} \
            {input.rna_genes} \
            {input.rna_barcodes} \
            {input.annot} \
            {input.gid} \
            {output.rna}
        
        echo "[INFO] Processing ATAC peaks to h5ad..."
        python workflow/scripts/dts/retina/process_peaks.py \
            {input.peaks_mtx} \
            {input.peak_names} \
            {input.peak_barcodes} \
            {input.annot} \
            {output.peaks}
        
        echo "[CLEANUP] Removing intermediate files..."
        rm -f {input.rna_mtx} \
              {input.rna_genes} \
              {input.rna_barcodes} \
              {input.peaks_mtx} \
              {input.peak_names} \
              {input.peak_barcodes}
        
        echo "[INFO] Retina download and processing complete!"
        """


rule annotate_retina:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        rna=rules.download_retina.output.rna,
        peaks=rules.download_retina.output.peaks,
        annot=rules.extract_retina_rds.output.annot,
    output: out='dts/retina/annotated.h5mu'
    shell:
        """
        python workflow/scripts/dts/retina/retina.py \
        -a {input.rna} \
        -b {input.peaks} \
        -c {input.annot} \
        -d {output.out}
        """
