rule compute_tfact_knocktf:
    input:
        grn='datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.{mdl}.grn.csv',
        rsc='/mnt/sds-hd/sd22b002/projects/GRETA/greta_resources/database/hg38/perturb/knocktf/',
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        temp(local('analysis/metrics/mech/tfact/knocktf/{dataset}.{case}.{pre}.{p2g}.{tfb}.{mdl}.scores.csv'))
    params:
        cats='config/prior_cats.json',
    shell:
        """
        python workflow/scripts/analysis/metrics/mech/compute_tfact.py \
        -i {input.grn} \
        -b {input.rsc} \
        -c {params.cats} \
        -o {output}
        """


rule compute_prtrb_knocktf:
    input:
        grn='datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.{mdl}.grn.csv',
        rsc='/mnt/sds-hd/sd22b002/projects/GRETA/greta_resources/database/hg38/perturb/knocktf/',
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        temp(local('analysis/metrics/mech/prtrb/knocktf/{dataset}.{case}.{pre}.{p2g}.{tfb}.{mdl}.scores.csv'))
    params:
        cats='config/prior_cats.json',
    resources:
        mem_mb=64000,
        runtime=1440,
    shell:
        """
        python workflow/scripts/analysis/metrics/mech/compute_prtrb.py \
        -i {input.grn} \
        -b {input.rsc} \
        -c {params.cats} \
        -o {output}
        """
