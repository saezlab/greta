rule download_reprofibro:
    output:
        tar=temp(local('datasets/reprofibro/RAW.tar')),
        barcodes=temp(local('datasets/reprofibro/barcode_map.tsv.gz')),
        genes=temp(local('datasets/reprofibro/genes.tsv.gz')),
        d1m_barcodes=temp(local('datasets/reprofibro/D1M.barcodes.tsv.gz')),
        d2m_barcodes=temp(local('datasets/reprofibro/D2M.barcodes.tsv.gz')),
        d1m_frags='datasets/reprofibro/D1M.frag.bed.gz',
        d2m_frags='datasets/reprofibro/D2M.frag.bed.gz',
        d1m_gex=temp(local('datasets/reprofibro/D1M.matrix.mtx.gz')),
        d2m_gex=temp(local('datasets/reprofibro/D2M.matrix.mtx.gz')),
        annot='datasets/reprofibro/annot.csv',
    params:
        tar=config['datasets']['reprofibro']['url']['tar'],
        barcodes=config['datasets']['reprofibro']['url']['barcodes'],
        genes=config['datasets']['reprofibro']['url']['genes'],
        annot=config['datasets']['reprofibro']['url']['annot'],
    shell:
        """
        data_path=$(dirname {output.tar})
        wget '{params.annot}' -O {output.annot}
        python workflow/scripts/datasets/reprofibro/prc_annot.py -a {output.annot}
        wget '{params.tar}' -O {output.tar}
        tar xvf {output.tar} -C $data_path
        mv $data_path/GSM7763381_D1M.barcodes.tsv.gz {output.d1m_barcodes}
        mv $data_path/GSM7763382_D2M.barcodes.tsv.gz {output.d2m_barcodes}
        mv $data_path/GSM7763381_D1M.frag.bed.gz {output.d1m_frags}
        mv $data_path/GSM7763382_D2M.frag.bed.gz {output.d2m_frags}
        mv $data_path/GSM7763381_D1M.matrix.mtx.gz {output.d1m_gex}
        mv $data_path/GSM7763382_D2M.matrix.mtx.gz {output.d2m_gex}
        wget '{params.barcodes}' -O {output.barcodes}
        wget '{params.genes}' -O {output.genes}
        """


rule callpeaks_reprofibro:
    threads: 32
    input:
        frags=['datasets/reprofibro/D1M.frag.bed.gz', 'datasets/reprofibro/D2M.frag.bed.gz'],
        annot='datasets/reprofibro/annot.csv',
    output:
        tmp=temp(directory(local('datasets/reprofibro/tmp'))),
        peaks=temp(local('datasets/reprofibro/peaks.h5ad'))
    resources:
        mem_mb=64000,
    shell:
        """
        python workflow/scripts/datasets/callpeaks.py \
        -f {input.frags} \
        -a {input.annot} \
        -t {output.tmp} \
        -o {output.peaks}
        """


rule annotate_reprofibro:
    input:
        path_matrix_d1m='datasets/reprofibro/D1M.matrix.mtx.gz',
        path_barcodes_d1m='datasets/reprofibro/D1M.barcodes.tsv.gz',
        path_matrix_d2m='datasets/reprofibro/D2M.matrix.mtx.gz',
        path_barcodes_d2m='datasets/reprofibro/D2M.barcodes.tsv.gz',
        path_gsym='datasets/reprofibro/genes.tsv.gz',
        path_peaks='datasets/reprofibro/peaks.h5ad',
        path_annot='datasets/reprofibro/annot.csv',
        path_barmap='datasets/reprofibro/barcode_map.tsv.gz',
        path_geneids='gdata/geneids/',
    output:
        'datasets/reprofibro/annotated.h5mu'
    params:
        organism=config['datasets']['reprofibro']['organism']
    shell:
        """
        python workflow/scripts/datasets/reprofibro/reprofibro.py \
        -a {input.path_matrix_d1m} \
        -b {input.path_barcodes_d1m} \
        -c {input.path_matrix_d2m} \
        -d {input.path_barcodes_d2m} \
        -e {input.path_gsym} \
        -f {input.path_peaks} \
        -g {input.path_annot} \
        -i {input.path_barmap} \
        -j {input.path_geneids} \
        -k {params.organism} \
        -l {output}
        """
