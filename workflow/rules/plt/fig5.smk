localrules: plt_fig5


rule plt_fig5:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: 'anl/dbs/stats.csv'
    output: 'plt/fig5/fig5.pdf'
    shell:
        """
        python workflow/scripts/plt/fig5/stats.py {input} {output}
        """
