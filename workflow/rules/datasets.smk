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
        organism='hg38',
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

# PBMC10K
rule download_pbmc10k:
    output:
        matrix=temp(local('datasets/pbmc10k/filtered_feature_bc_matrix.h5')),
        peak_annot=temp(local('datasets/pbmc10k/atac_peak_annotation.tsv')),
        atac_frags=temp(local('datasets/pbmc10k/atac_fragments.tsv.gz')),
        atac_index=temp(local('datasets/pbmc10k/atac_fragments.tsv.gz.tbi')),
    params:
        matrix=config['datasets']['pbmc10k']['url']['matrix'],
        peak_annot=config['datasets']['pbmc10k']['url']['peak_annot'],
        atac_frags=config['datasets']['pbmc10k']['url']['atac_frags'],
        atac_index=config['datasets']['pbmc10k']['url']['atac_index']
    shell:
        """
        wget '{params.matrix}' -O {output.matrix}
        wget '{params.peak_annot}' -O {output.peak_annot}
        wget '{params.atac_frags}' -O {output.atac_frags}
        wget '{params.atac_index}' -O {output.atac_index}
        """
