rule mech_tfact:
    input:
        grn=lambda wildcards: rules.grn_run.output.out.format(**wildcards),
        rsc='/mnt/sds-hd/sd22b002/projects/GRETA/greta_resources/database/hg38/perturb/{resource}/',
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        out='analysis/metrics/mech/tfact/{resource}/{dataset}.{case}/{pre}.{p2g}.{tfb}.{mdl}.scores.csv'
    params:
        cats='config/prior_cats.json',
    shell:
        """
        python workflow/scripts/analysis/metrics/mech/compute_tfact.py \
        -i {input.grn} \
        -b {input.rsc} \
        -c {params.cats} \
        -o {output.out}
        """


rule mech_prtrb:
    input:
        grn=lambda wildcards: rules.grn_run.output.out.format(**wildcards),
        rsc='/mnt/sds-hd/sd22b002/projects/GRETA/greta_resources/database/hg38/perturb/{resource}/',
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        out='analysis/metrics/mech/prtrb/{resource}/{dataset}.{case}/{pre}.{p2g}.{tfb}.{mdl}.scores.csv'
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
        -o {output.out}
        """
