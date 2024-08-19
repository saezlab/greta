# heartatlas
rule download_heartatlas:
    output:
        tar=local('datasets/heartatlas/fragments.tar'),
        frag=local(expand('datasets/heartatlas/{batchID}_atac_fragments.tsv.gz', batchID=config['datasets']['heartatlas']['batchIDs'])),
        anndata=local('datasets/heartatlas/multiome_raw.h5ad')
    params:
        tar=config['datasets']['heartatlas']['url']['tar'],
        anndata=config['datasets']['heartatlas']['url']['anndata']
    shell:
        """
        data_path=$(dirname {output.tar})
        wget '{params.tar}' -O '{output.tar}'
        tar -xvf {output.tar} -C $data_path
        rm $data_path/*.tbi
        wget '{params.anndata}' -O '{output.anndata}'
        """

rule prcannot_heartatlas:
    input:
        h5ad='datasets/heartatlas/multiome_raw.h5ad'
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        annot='datasets/heartatlas/annot.csv'
    shell:
        """
        python workflow/scripts/datasets/heartatlas/heart_annot.py \
        -h {input.h5ad} \
        -a {output.annot}
        """

rule callpeaks_heartatlas:
    input:
        frags=expand('datasets/heartatlas/{batchID}_atac_fragments.tsv.gz', batchID=config['datasets']['heartatlas']['batchIDs']),
        annot='datasets/heartatlas/annot.csv',
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        tmp=directory(local('datasets/heartatlas/tmp_peaks')),
        peaks=local('datasets/heartatlas/peaks.h5ad')
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