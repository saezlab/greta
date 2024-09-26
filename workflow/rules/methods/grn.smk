rule grn_run:
    input:
        mdl=lambda wildcards: map_rules('mdl', wildcards.mdl),
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        out='datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.{mdl}.grn.csv'
    resources:
        mem_mb=32000,
    shell:
        """
        python workflow/scripts/methods/grn.py \
        -i {input.mdl} \
        -o {output.out}
        """


rule mdl_collectri:
    input:
        mdata=rules.extract_case.output.mdata,
        grn='/mnt/sds-hd/sd22b002/projects/GRETA/greta_resources/database/hg38/grn/collectri.csv',
        proms='/mnt/sds-hd/sd22b002/projects/GRETA/greta_resources/database/hg38/cre/promoters/promoters.bed',
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        out='datasets/{dataset}/cases/{case}/runs/collectri.collectri.collectri.collectri.mdl.csv'
    resources:
        mem_mb=32000
    shell:
        """
        python workflow/scripts/methods/prc_prior_grn.py \
        -g {input.grn} \
        -d {input.mdata} \
        -p {input.proms} \
        -o {output.out}
        """


rule mdl_dorothea:
    input:
        mdata=rules.extract_case.output.mdata,
        grn='/mnt/sds-hd/sd22b002/projects/GRETA/greta_resources/database/hg38/grn/dorothea.csv',
        proms='/mnt/sds-hd/sd22b002/projects/GRETA/greta_resources/database/hg38/cre/promoters/promoters.bed',
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        out='datasets/{dataset}/cases/{case}/runs/dorothea.dorothea.dorothea.dorothea.mdl.csv'
    resources:
        mem_mb=32000
    shell:
        """
        python workflow/scripts/methods/prc_prior_grn.py \
        -g {input.grn} \
        -d {input.mdata} \
        -p {input.proms} \
        -o {output.out}
        """
