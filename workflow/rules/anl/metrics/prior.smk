localrules: prior_tfm, prior_cre, prior_eqtl


rule prior_tfm:
    input:
        grn=lambda wildcards: rules.grn_run.output.out.format(**wildcards),
        db='dbs/hg38/tfm/{db}/{db}.tsv',
    params:
        cats='config/prior_cats.json',
    output:
        out='anl/metrics/prior/tfm/{db}/{dat}.{case}/{pre}.{p2g}.{tfb}.{mdl}.scores.csv'
    shell:
        """
        python workflow/scripts/anl/metrics/prior/compute_tfm.py \
        -a {input.grn} \
        -b {input.db} \
        -c {params.cats} \
        -f {output.out}
        """


rule prior_tfbind:
    input:
        grn=lambda wildcards: rules.grn_run.output.out.format(**wildcards),
        db='dbs/hg38/tfb/{db}/{db}.bed',
    params:
        cats='config/prior_cats.json',
    output:
        out='anl/metrics/prior/tfb/{db}/{dat}.{case}/{pre}.{p2g}.{tfb}.{mdl}.scores.csv'
    params:
        grp='source',
    shell:
        """
        python workflow/scripts/anl/metrics/prior/compute_gnm.py \
        -a {input.grn} \
        -b {input.db} \
        -c {params.cats} \
        -d {params.grp} \
        -f {output}
        """


rule prior_cre:
    input:
        grn=lambda wildcards: rules.grn_run.output.out.format(**wildcards),
        db='dbs/hg38/cre/{db}/{db}.bed',
    params:
        cats='config/prior_cats.json',
    output:
        out='anl/metrics/prior/cre/{db}/{dat}.{case}/{pre}.{p2g}.{tfb}.{mdl}.scores.csv'
    shell:
        """
        python workflow/scripts/anl/metrics/prior/compute_gnm.py \
        -a {input.grn} \
        -b {input.resource} \
        -c {params.cats} \
        -f {output}
        """


rule prior_eqtl:
    input:
        grn=lambda wildcards: rules.grn_run.output.out.format(**wildcards),
        resource='dbs/hg38/eqtl/{db}/{db}.bed',
    params:
        cats='config/prior_cats.json',
    output:
        out='anl/metrics/prior/eqtl/{db}/{dat}.{case}/{pre}.{p2g}.{tfb}.{mdl}.scores.csv'
    params:
        grp='target',
    shell:
        """
        python workflow/scripts/anl/metrics/prior/compute_gnm.py \
        -a {input.grn} \
        -b {input.resource} \
        -c {params.cats} \
        -d {params.grp} \
        -f {output}
        """
