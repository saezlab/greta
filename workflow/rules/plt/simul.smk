localrules: fig_simul

rule fig_simul:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: rules.topo_simul.output
    output: 'plt/simul/fig.pdf'
    shell:
        """
        python workflow/scripts/plt/simul/simul.py {input} {output}
        """
