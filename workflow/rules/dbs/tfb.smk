localrules: tfb_m_chipatlas, tfb_t_chipatlas, tfb_chipatlas, tfb_m_remap2022, tfb_t_remap2022, tfb_remap2022, tfb_t_unibind, tfb_unibind


checkpoint tfb_m_chipatlas:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: rules.lambert.output[0]
    output: 'dbs/hg38/tfb/chipatlas/meta.tsv'
    params:
        mta=config['dbs']['hg38']['tfb']['chipatlas']['meta'],
    shell:
        """
        wget --no-verbose '{params.mta}' -O {output} && \
        python workflow/scripts/dbs/tfb/chipatlas_meta.py {output} {input}
        """


rule tfb_t_chipatlas:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: 'dbs/hg38/tfb/chipatlas/meta.tsv'
    output: 'dbs/hg38/tfb/chipatlas/raw/{chipatlas_tf}.bed'
    params:
        url=config['dbs']['hg38']['tfb']['chipatlas']['url'],
        max_psize=config['dbs']['hg38']['tfb']['max_psize']
    shell:
        """
        if wget --spider '{params.url}' 2>/dev/null; then
            wget --no-verbose --tries=2 '{params.url}' -O - | \
            python workflow/scripts/dbs/tfb/chipatlas_tf.py {output} {input} {params.max_psize} | \
            bedtools merge -c 4,5 -o distinct,distinct > {output}
        else
            touch {output}
        fi
        """


def chipatlas_aggr(wildcards):
    import pandas as pd
    meta_file = checkpoints.tfb_m_chipatlas.get().output[0]
    tfs = sorted(pd.read_csv(meta_file, sep='\t', header=None)[1].unique())
    return expand("dbs/hg38/tfb/chipatlas/raw/{chipatlas_tf}.bed", chipatlas_tf=tfs)


rule tfb_chipatlas:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: chipatlas_aggr
    output: 'dbs/hg38/tfb/chipatlas/chipatlas.bed'
    shell:
        """
        cat dbs/hg38/tfb/chipatlas/raw/*.bed |
        python workflow/scripts/dbs/tfb/aggregate.py > {output}
        """


rule tfb_m_remap2022:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: rules.ont_bto.output[0]
    output: 'dbs/hg38/tfb/remap2022/meta.tsv'
    params:
        mta=config['dbs']['hg38']['tfb']['remap2022']['meta'],
    shell:
        """
        wget --no-verbose '{params.mta}' -O - | \
        python workflow/scripts/dbs/tfb/remap2022_meta.py {input} {output}
        """


checkpoint tfb_r_remap2022:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        tfs=rules.lambert.output[0],
        mta=rules.tfb_m_remap2022.output
    output: directory('dbs/hg38/tfb/remap2022/raw')
    params:
        url=config['dbs']['hg38']['tfb']['remap2022']['url'],
        max_psize=config['dbs']['hg38']['tfb']['max_psize']
    shell:
        """
        mkdir -p {output} && \
        wget --no-verbose '{params.url}' -O {output}.tmp && \
        zcat {output}.tmp | \
        python workflow/scripts/dbs/tfb/remap2022_raw.py \
        {input.tfs} \
        {input.mta} \
        {params.max_psize} \
        {output} && \
        rm {output}.tmp
        """


rule tfb_t_remap2022:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: 'dbs/hg38/tfb/remap2022/raw/{remap2022_tf}.bed'
    output: 'dbs/hg38/tfb/remap2022/raw/{remap2022_tf}.bed'
    shell:
        """
        bedtools merge -i {output} -c 4,5 -o distinct,distinct > {output}.tmp && \
        mv {output}.tmp {output}
        """


def remap2022_aggr(wildcards):
    checkpoints.tfb_r_remap2022.get()
    tfs = glob_wildcards(f"dbs/hg38/tfb/remap2022/raw/{{tf}}.bed").tf
    return expand(f"dbs/hg38/tfb/remap2022/raw/{tf}.bed", tf=tfs)


rule tfb_remap2022:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: remap2022_aggr
    output: 'dbs/hg38/tfb/remap2022/remap2022.bed'
    shell:
        """
        cat dbs/hg38/tfb/remap2022/raw/*.bed |
        python workflow/scripts/dbs/tfb/aggregate.py > {output}
        """


checkpoint tfb_r_unibind:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: rules.lambert.output,
    output: directory('dbs/hg38/tfb/unibind/raw')
    params:
        url=config['dbs']['hg38']['tfb']['unibind']['url'],
        max_psize=config['dbs']['hg38']['tfb']['max_psize']
    shell:
        """
        mkdir -p {output} && \
        wget --no-verbose '{params.url}' -O {output}.tmp && \
        zcat {output}.tmp | \
        python workflow/scripts/dbs/tfb/unibind_raw.py \
        {input} \
        {params.max_psize} \
        {output} && \
        rm {output}.tmp
        """


rule tfb_t_unibind:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: 'dbs/hg38/tfb/unibind/raw/{unibind_tf}.bed'
    output: 'dbs/hg38/tfb/unibind/raw/{unibind_tf}.bed'
    shell:
        """
        bedtools merge -i {output} -c 4,5 -o distinct,distinct > {output}.tmp && \
        mv {output}.tmp {output}
        """


def unibind_aggr(wildcards):
    checkpoints.tfb_r_unibind.get()
    tfs = glob_wildcards(f"dbs/hg38/tfb/unibind/raw/{{tf}}.bed").tf
    return expand(f"dbs/hg38/tfb/unibind/raw/{tf}.bed", tf=tfs)


rule tfb_unibind:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: unibind_aggr
    output: 'dbs/hg38/tfb/unibind/unibind.bed'
    shell:
        """
        cat dbs/hg38/tfb/unibind/raw/*.bed |
        python workflow/scripts/dbs/tfb/aggregate.py > {output}
        """
