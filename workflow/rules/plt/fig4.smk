localrules: plt_fig4


rule plt_fig4:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        sims='anl/topo/pbmc10k.all.sims_mult.csv',
        stat='anl/topo/pbmc10k.all.stats_mult.csv',
    output: 'plt/fig4/fig4.pdf'
    shell:
        """
        python workflow/scripts/plt/fig4/sims.py {input.sims} {input.stat} {output}
        """
