localrules: prior_tfm


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



rule all:
    input:
        "results/f01_results.csv"

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
        cont=temp(local("results/{grn}_contingency.csv")),
        db="data/{database}.csv"
    output:
        out=temp(local("results/{grn}_{database}_f01.csv"))
    shell:
        """
        python scripts/f01.py -c {input.cont} -d {input.db} -o {output.out}
        """

rule merge_results:
    input:
        temp(local(expand("results/{grn}_{database}_f01.csv",
               grn=["granie", "figr", "pando", "random", "collectri"],
               database=["intact", "pubmed"])))
    output:
        "results/f01_results.csv"
    shell:
        """
        python -c "import pandas as pd; import sys; \
        files = sys.argv[1:]; \
        df = pd.concat([pd.read_csv(file).assign(Source=file.split('/')[-1].replace('.csv', '')) for file in files]); \
        df.to_csv('{output}', index=False);" {input}
        """











