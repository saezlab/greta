rule get_grn:
    input:
        'datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.{mdl}.mdl.csv'
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        'datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.{mdl}.grn.csv'
    resources:
        mem_mb=16000,
    shell:
        """
        python workflow/scripts/methods/grn.py \
        -i {input} \
        -o {output}
        """

rule prc_prior_grn:
    input:
        data='datasets/{dataset}/cases/{case}/mdata.h5mu',
        grn='/mnt/sds-hd/sd22b002/projects/GRETA/greta_resources/database/hg38/grn/{grn_name}.csv',
        proms='/mnt/sds-hd/sd22b002/projects/GRETA/greta_resources/database/hg38/cre/promoters/promoters.bed',
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        'datasets/{dataset}/cases/{case}/runs/{grn_name}.grn.csv'
    resources:
        mem_mb=16000
    shell:
        """
        python workflow/scripts/methods/prc_prior_grn.py \
        -g {input.grn} \
        -d {input.data} \
        -p {input.proms} \
        -o {output}
        """
