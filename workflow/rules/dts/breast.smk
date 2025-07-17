localrules: prcannot_breast

rule download_anndata_breast:
    threads: 1
    output:
        adata=temp(local('dts/breast/multiome_raw.h5ad')),
    params:
        adata=config['dts']['breast']['url']['anndata']
    shell:
        """
        wget --no-verbose '{params.adata}' -O '{output.adata}'
        """


rule prcannot_breast:
    threads: 1
    singularity:
        'workflow/envs/gretabench.sif'
    input: rules.download_anndata_breast.output.adata,
    output:
        annot=temp(local('dts/breast/annot.csv'))
    shell:
        """
        python workflow/scripts/dts/breast/breast_annot.py \
        -i {input} \
        -o {output.annot}
        """