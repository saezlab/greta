localrules: genom_cre

rule genom_tfb:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        grn=lambda wildcards: rules.grn_run.output.out.format(**wildcards),
        db='dbs/{org}/tfb/{db}/{db}.bed',
    output:
        out='anl/metrics/genom/tfb/{db}/{org}.{dat}.{case}/{pre}.{p2g}.{tfb}.{mdl}.scores.csv'
    params:
        grp='source',
    shell:
        """
        python workflow/scripts/anl/metrics/genom/gnm.py \
        -a {input.grn} \
        -b {input.db} \
        -d {params.grp} \
        -f {output}
        """


rule genom_cre:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        grn=lambda wildcards: rules.grn_run.output.out.format(**wildcards),
        db='dbs/{org}/cre/{db}/{db}.bed',
    output:
        out='anl/metrics/genom/cre/{db}/{org}.{dat}.{case}/{pre}.{p2g}.{tfb}.{mdl}.scores.csv'
    shell:
        """
        python workflow/scripts/anl/metrics/genom/gnm.py \
        -a {input.grn} \
        -b {input.db} \
        -f {output}
        """


rule genom_c2g:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        grn=lambda wildcards: rules.grn_run.output.out.format(**wildcards),
        resource='dbs/{org}/c2g/{db}/{db}.bed',
    output:
        out='anl/metrics/genom/c2g/{db}/{org}.{dat}.{case}/{pre}.{p2g}.{tfb}.{mdl}.scores.csv'
    params:
        grp='target',
    shell:
        """
        python workflow/scripts/anl/metrics/genom/gnm.py \
        -a {input.grn} \
        -b {input.resource} \
        -d {params.grp} \
        -f {output}
        """