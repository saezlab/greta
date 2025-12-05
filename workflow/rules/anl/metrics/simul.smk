localrules: simul_simul

rule simul_simul:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: rules.topo_simul.output
    output: 'anl/metrics/simul/simul.scores.csv'
    shell:
        """
        python workflow/scripts/anl/metrics/simul/simul.py {output}
        """
