localrules: prior_tfm, prior_tfp, prior_grn


rule prior_tfm:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        grn=lambda wildcards: rules.grn_run.output.out.format(**wildcards),
        db='dbs/{org}/tfm/{db}/{db}.tsv.gz',
    output:
        out='anl/metrics/prior/tfm/{db}/{org}.{dat}.{case}/{pre}.{p2g}.{tfb}.{mdl}.scores.csv'
    shell:
        """
        python workflow/scripts/anl/metrics/prior/tfm.py \
        -a {input.grn} \
        -b {input.db} \
        -f {output.out}
        """


rule prior_tfp:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        grn=lambda wildcards: rules.grn_run.output.out.format(**wildcards),
        db='dbs/{org}/tfp/{db}/{db}.tsv.gz',
    output:
        out='anl/metrics/prior/tfp/{db}/{org}.{dat}.{case}/{pre}.{p2g}.{tfb}.{mdl}.scores.csv'
    params:
        thr_p=0.01,
    shell:
        """
        python workflow/scripts/anl/metrics/prior/tfp.py \
        {input.grn} {input.db} {params.thr_p} {output.out}
        """


rule prior_grn:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        grn=lambda wildcards: rules.grn_run.output.out.format(**wildcards),
        db='dbs/{org}/gst/{db}.csv.gz',
    output:
        out='anl/metrics/prior/grn/{db}/{org}.{dat}.{case}/{pre}.{p2g}.{tfb}.{mdl}.scores.csv'
    shell:
        """
        python workflow/scripts/anl/metrics/prior/grn.py \
        -a {input.grn} \
        -b {input.db} \
        -f {output.out}
        """
