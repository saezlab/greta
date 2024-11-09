localrules: tfb_m_chipatlas, tfb_t_chipatlas, tfb_chipatlas

# chipatlas  remap2022  unibind  union

checkpoint tfb_m_chipatlas:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: rules.lambert.output[0]
    output: 'dbs/hg38/tfb/chipatlas/raw/meta.tsv'
    params:
        mta=config['dbs']['hg38']['tfb']['chipatlas']['meta'],
        #url=config['dbs']['hg38']['tfb']['chipatlas']['url']
    shell:
        """
        wget --no-verbose '{params.mta}' -O {output} && \
        python workflow/scripts/dbs/tfb/chipatlas_meta.py {output} {input}
        """


rule tfb_t_chipatlas:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: 'dbs/hg38/tfb/chipatlas/raw/meta.tsv'
    output: 'dbs/hg38/tfb/chipatlas/raw/{chipatlas_tf}.bed'
    params:
        url=config['dbs']['hg38']['tfb']['chipatlas']['url'],
        max_psize=config['dbs']['hg38']['tfb']['max_psize']
    resources:
        mem_mb=2000
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
        python workflow/scripts/dbs/tfb/chipatlas.py > {output}
        """



