localrules: prior_tfm, prior_cre, prior_c2g


rule prior_tfm:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        grn=lambda wildcards: rules.grn_run.output.out.format(**wildcards),
        db='dbs/hg38/tfm/{db}/{db}.tsv',
    params:
        cats='config/prior_cats.json',
    output:
        out='anl/metrics/prior/tfm/{db}/{dat}.{case}/{pre}.{p2g}.{tfb}.{mdl}.scores.csv'
    shell:
        """
        python workflow/scripts/anl/metrics/prior/tfm.py \
        -a {input.grn} \
        -b {input.db} \
        -c {params.cats} \
        -f {output.out}
        """


rule prior_tfb:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        grn=lambda wildcards: rules.grn_run.output.out.format(**wildcards),
        db='dbs/hg38/tfb/{db}/{db}.bed',
    output:
        out='anl/metrics/prior/tfb/{db}/{dat}.{case}/{pre}.{p2g}.{tfb}.{mdl}.scores.csv'
    params:
        cats='config/prior_cats.json',
        grp='source',
    shell:
        """
        python workflow/scripts/anl/metrics/prior/gnm.py \
        -a {input.grn} \
        -b {input.db} \
        -c {params.cats} \
        -d {params.grp} \
        -f {output}
        """


rule prior_cre:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        grn=lambda wildcards: rules.grn_run.output.out.format(**wildcards),
        db='dbs/hg38/cre/{db}/{db}.bed',
    params:
        cats='config/prior_cats.json',
    output:
        out='anl/metrics/prior/cre/{db}/{dat}.{case}/{pre}.{p2g}.{tfb}.{mdl}.scores.csv'
    shell:
        """
        python workflow/scripts/anl/metrics/prior/gnm.py \
        -a {input.grn} \
        -b {input.db} \
        -c {params.cats} \
        -f {output}
        """


rule prior_c2g:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        grn=lambda wildcards: rules.grn_run.output.out.format(**wildcards),
        resource='dbs/hg38/c2g/{db}/{db}.bed',
    output:
        out='anl/metrics/prior/c2g/{db}/{dat}.{case}/{pre}.{p2g}.{tfb}.{mdl}.scores.csv'
    params:
        cats='config/prior_cats.json',
        grp='target',
    shell:
        """
        python workflow/scripts/anl/metrics/prior/gnm.py \
        -a {input.grn} \
        -b {input.resource} \
        -c {params.cats} \
        -d {params.grp} \
        -f {output}
        """
