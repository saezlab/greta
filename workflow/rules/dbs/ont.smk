localrules: ont_bto


checkpoint ont_bto:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    output: 'dbs/ont/bto.tsv'
    params:
        url=config['dbs']['ont']['bto'],
    shell:
        """
        wget --no-verbose '{params.url}' -O - | \
        python workflow/scripts/dbs/ont/bto.py {output}
        """