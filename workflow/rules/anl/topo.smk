rule topo_mult:
    threads: 32
    singularity: 'workflow/envs/gretabench.sif'
    input:
        lambda w: make_combs_rules(w=w, mthds=mthds, baselines=baselines, rule_name='grn_run')
    output:
        stats='anl/topo/{dat}.{case}.stats_mult.csv',
        sims='anl/topo/{dat}.{case}.sims_mult.csv',
    shell:
        """
        python workflow/scripts/anl/topo/run_pair_sim.py \
        -p {input} \
        -t {output.stats} \
        -s {output.sims}
        """
