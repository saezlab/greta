localrules: prior_tfm, prior_tfp


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


rule prior_tfp:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        grn=lambda wildcards: rules.grn_run.output.out.format(**wildcards),
        db='dbs/hg38/tfp/{db}/{db}.tsv',
    output:
        out='anl/metrics/prior/tfp/{db}/{dat}.{case}/{pre}.{p2g}.{tfb}.{mdl}.scores.csv'
    params:
        thr_p=0.01,
    shell:
        """
        python workflow/scripts/anl/metrics/prior/tfp.py \
        {input.grn} {input.db} {params.thr_p} {output.out}
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


rule calculate_contingency:
    input:
        tf="data/lambert.csv",
        grn="data/{grn}.csv"
    output:
        tmp_cont=temp(local("results/{grn}_contingency.csv"))
    shell:
        """
        python scripts/cont.py -t {input.tf} -m {input.grn} -o {output.tmp_cont}
        """

rule calculate_f01:
    input:
        cont="results/{grn}_contingency.csv",
        db="data/{database}.csv"
    output:
        out=temp(local("results/{grn}_{database}_f01.csv"))
    shell:
        """
        python scripts/f01.py -c {input.cont} -d {input.db} -o {output.out}
        """


