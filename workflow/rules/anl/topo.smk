localrules: topo_inter


rule topo_mult:
    threads: 4
    singularity: 'workflow/envs/gretabench.sif'
    input:
        lambda w: make_combs_rules(w=w, mthds=mthds, baselines=baselines, rule_name='grn_run')
    output:
        stats='anl/topo/{dat}.{case}.stats_mult.csv',
        sims='anl/topo/{dat}.{case}.sims_mult.csv',
    resources:
        mem_mb=128000
    shell:
        """
        python workflow/scripts/anl/topo/run_pair_sim.py \
        -t {output.stats} \
        -s {output.sims}
        """


rule topo_inter:
    threads: 1
    input:
        lambda w: make_combs_rules(w=w, mthds=mthds, baselines=baselines, rule_name='grn_run')
    output: 'anl/topo/{dat}.{case}.inter.csv',
    params: min_prop=config['topo_min_prop']
    shell:
        """
        python workflow/scripts/anl/topo/inter.py \
        -g {input} \
        -b {baselines} \
        -p {params.min_prop} \
        -o {output}
        """
