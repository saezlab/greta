rule extract_case:
    input:
        'datasets/{dataset}/annotated.h5mu'
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        mdata='datasets/{dataset}/cases/{case}/mdata.h5mu',
    params:
        celltypes=lambda w: config['datasets'][w.dataset]['cases'][w.case]['celltypes'],
        n_hvg=lambda w: config['datasets'][w.dataset]['cases'][w.case]['n_hvg'],
        n_hvr=lambda w: config['datasets'][w.dataset]['cases'][w.case]['n_hvr'],
    shell:
        """
        python workflow/scripts/datasets/extract_case.py \
        -i '{input}' \
        -c '{params.celltypes}' \
        -g '{params.n_hvg}' \
        -r '{params.n_hvr}' \
        -o '{output.mdata}'
        """

rule download_geneids:
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        hg='gdata/geneids/hg38.csv',
        mm='gdata/geneids/mm10.csv',
        dr=directory('gdata/geneids')
    shell:
        """
        Rscript workflow/scripts/datasets/download_geneids.R \
        {output.hg} \
        {output.mm}
        """

# NEURIPS2021
rule download_neurips2021:
    output:
        temp(local('datasets/neurips2021/original.h5ad'))
    params:
        url=config['datasets']['neurips2021']['url']
    shell:
        """
        wget '{params.url}' -O '{output}.gz'
        gzip -d {output}.gz
        """

rule annotate_neurips2021:
    input:
        d='datasets/neurips2021/original.h5ad',
        g='gdata/geneids',
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        'datasets/neurips2021/annotated.h5mu'
    params:
        organism=config['datasets']['neurips2021']['organism'],
    resources:
        mem_mb=48000,
    shell:
        """
        python workflow/scripts/datasets/neurips2021.py \
        -i {input.d} \
        -r {params.organism} \
        -g {input.g} \
        -o {output}
        """


# reprofibro
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
    input:
        frags=['datasets/reprofibro/D1M.frag.bed.gz', 'datasets/reprofibro/D2M.frag.bed.gz'],
        annot='datasets/reprofibro/annot.csv',
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        tmp=temp(directory(local('datasets/reprofibro/tmp'))),
        peaks=temp(local('datasets/reprofibro/peaks.h5ad'))
    resources:
        mem_mb=64000,
    threads: 16
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
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        'datasets/reprofibro/annotated.h5mu'
    params:
        organism=config['datasets']['reprofibro']['organism'],
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


# PBMC10K
rule download_pbmc10k:
    output:
        atac_frags='datasets/pbmc10k/smpl.frags.tsv.gz',
    params:
        matrix=config['datasets']['pbmc10k']['url']['matrix'],
        atac_frags=config['datasets']['pbmc10k']['url']['atac_frags'],
    shell:
        """
        wget '{params.atac_frags}' -O {output.atac_frags}
        """

rule prcannot_pbmc10k:
    output:
        tmp=temp(directory(local('datasets/pbmc10k/tmp'))),
        annot=temp(local('datasets/pbmc10k/annot.csv')),
    shell:
        """
        python workflow/scripts/datasets/pbmc10k/prc_annot.py \
        -t {output.tmp} \
        -a {output.annot}
        """

rule callpeaks_pbmc10k:
    input:
        frags='datasets/pbmc10k/smpl.frags.tsv.gz',
        annot='datasets/pbmc10k/annot.csv',
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        tmp=temp(directory(local('datasets/pbmc10k/tmp_peaks'))),
        peaks=temp(local('datasets/pbmc10k/peaks.h5ad'))
    resources:
        mem_mb=64000,
    threads: 16
    shell:
        """
        python workflow/scripts/datasets/callpeaks.py \
        -f {input.frags} \
        -a {input.annot} \
        -t {output.tmp} \
        -o {output.peaks}
        """

rule annotate_pbmc10k:
    input:
        annot='datasets/pbmc10k/annot.csv',
        g='gdata/geneids',
        peaks='datasets/pbmc10k/peaks.h5ad',
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        tmp=temp(directory(local('datasets/pbmc10k/tmp_annot'))),
        out='datasets/pbmc10k/annotated.h5mu'
    params:
        organism=config['datasets']['pbmc10k']['organism'],
    resources:
        mem_mb=32000,
    shell:
        """
        python workflow/scripts/datasets/pbmc10k/pbmc10k.py \
        -a {output.tmp} \
        -b {input.annot} \
        -c {input.g} \
        -d {params.organism} \
        -e {input.peaks} \
        -f {output.out}
        """
