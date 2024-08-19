rule compute_tfact_knocktf:
    input:
        'datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.{mdl}.grn.csv'
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        temp(local('analysis/metrics/mech/tfact/knocktf/{dataset}.{case}.{pre}.{p2g}.{tfb}.{mdl}.scores.csv'))
    params:
        knocktf='/mnt/sds-hd/sd22b002/projects/GRETA/greta_resources/database/hg38/perturb/oldknocktf/',
        cats=['Blood', 'Haematopoietic_and_lymphoid_tissue', 'Haematopoietic_and_lymphoid_tissue_Blood']
    shell:
        """
        python workflow/scripts/analysis/metrics/mech/compute_tfact.py \
        -i {input} \
        -b {params.knocktf} \
        -c {params.cats} \
        -o {output}
        """


rule compute_prtrb_knocktf:
    input:
        'datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.{mdl}.grn.csv'
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        temp(local('analysis/metrics/mech/prtrb/knocktf/{dataset}.{case}.{pre}.{p2g}.{tfb}.{mdl}.scores.csv'))
    params:
        knocktf='/mnt/sds-hd/sd22b002/projects/GRETA/greta_resources/database/hg38/perturb/oldknocktf/',
        cats=['Blood', 'Haematopoietic_and_lymphoid_tissue', 'Haematopoietic_and_lymphoid_tissue_Blood']
    resources:
        mem_mb=64000,
        runtime=360,
    shell:
        """
        python workflow/scripts/analysis/metrics/mech/compute_prtrb.py \
        -i {input} \
        -b {params.knocktf} \
        -c {params.cats} \
        -o {output}
        """
