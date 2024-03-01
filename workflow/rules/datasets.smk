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

# NEURIPS2021
rule download_neurips2021:
    output:
        'datasets/neurips2021/original.h5ad.gz'
    params:
        url=config['datasets']['neurips2021']['url']
    shell:
        "wget '{params.url}' -O '{output}'"

rule decompress_neurips2021:
    input:
        'datasets/neurips2021/original.h5ad.gz'
    output:
        'datasets/neurips2021/original.h5ad'
    shell:
        'gzip -d {input}'
        
rule annotate_neurips2021:
    resources:
        mem_mb=48000,
    input: 'datasets/neurips2021/original.h5ad'
    output:
        mdata='datasets/neurips2021/annotated.h5mu'
    singularity:
        'workflow/envs/gretabench.sif'
    shell:
        'python workflow/scripts/datasets/neurips2021.py -i {input} -o {output.mdata}'

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
