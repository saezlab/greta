localrules: fig_dbs


rule fig_dbs:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        sts='anl/dbs/stats.csv',
        ovc='anl/dbs/ocoef.csv',
    output: 'plt/dbs/fig.pdf'
    shell:
        """
        python workflow/scripts/plt/dbs/stats.py {input.sts} {input.ovc} {output}
        """
