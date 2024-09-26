localrules: prior_tfm, prior_cre, prior_eqtl


rule prior_tfm:
    input:
        grn=lambda wildcards: rules.grn_run.output.out.format(**wildcards),
        resource='/mnt/sds-hd/sd22b002/projects/GRETA/greta_resources/database/hg38/tfm/{resource}/{resource}.csv',
    params:
        cats='config/prior_cats.json',
    output:
        out='analysis/metrics/prior/tfm/{resource}/{dataset}.{case}/{pre}.{p2g}.{tfb}.{mdl}.scores.csv'
    shell:
        """
        python workflow/scripts/analysis/metrics/prior/compute_tfm.py \
        -a {input.grn} \
        -b {input.resource} \
        -c {params.cats} \
        -f {output.out}
        """


rule prior_tfbind:
    input:
        grn=lambda wildcards: rules.grn_run.output.out.format(**wildcards),
        resource='/mnt/sds-hd/sd22b002/projects/GRETA/greta_resources/database/hg38/tfbind/{resource}/{resource}.bed',
    params:
        cats='config/prior_cats.json',
    output:
        out='analysis/metrics/prior/tfbind/{resource}/{dataset}.{case}/{pre}.{p2g}.{tfb}.{mdl}.scores.csv'
    params:
        grp='source',
    shell:
        """
        python workflow/scripts/analysis/metrics/prior/compute_gnm.py \
        -a {input.grn} \
        -b {input.resource} \
        -c {params.cats} \
        -d {params.grp} \
        -f {output}
        """


rule prior_cre:
    input:
        grn=lambda wildcards: rules.grn_run.output.out.format(**wildcards),
        resource='/mnt/sds-hd/sd22b002/projects/GRETA/greta_resources/database/hg38/cre/{resource}/{resource}.bed',
    params:
        cats='config/prior_cats.json',
    output:
        out='analysis/metrics/prior/cre/{resource}/{dataset}.{case}/{pre}.{p2g}.{tfb}.{mdl}.scores.csv'
    shell:
        """
        python workflow/scripts/analysis/metrics/prior/compute_gnm.py \
        -a {input.grn} \
        -b {input.resource} \
        -c {params.cats} \
        -f {output}
        """


rule prior_eqtl:
    input:
        grn=lambda wildcards: rules.grn_run.output.out.format(**wildcards),
        resource='/mnt/sds-hd/sd22b002/projects/GRETA/greta_resources/database/hg38/eqtl/{resource}/{resource}.bed',
    params:
        cats='config/prior_cats.json',
    output:
        out='analysis/metrics/prior/eqtl/{resource}/{dataset}.{case}/{pre}.{p2g}.{tfb}.{mdl}.scores.csv'
    params:
        grp='target',
    shell:
        """
        python workflow/scripts/analysis/metrics/prior/compute_gnm.py \
        -a {input.grn} \
        -b {input.resource} \
        -c {params.cats} \
        -d {params.grp} \
        -f {output}
        """
