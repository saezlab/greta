rule topo_mult:
    threads: 32
    input:
        lambda w: make_combs_rules(w=w, mthds=mthds, baselines=baselines, rule_name='grn_run')
    output:
        stats='analysis/topo/{dataset}.{case}.stats_mult.csv',
        sims='analysis/topo/{dataset}.{case}.sims_mult.csv',
    shell:
        """
        python workflow/scripts/analysis/topo/run_pair_sim.py \
        -p {input} \
        -t {output.stats} \
        -s {output.sims}
        """

rule topo_mult_pitu:
    threads: 32
    input:
        lambda w: make_combs_rules(w=w, mthds=mthds, baselines=baselines, rule_name='grn_run')
    output:
        stats='analysis/topo/pitupair.all.stats_mult.csv',
        sims='analysis/topo/pitupair.all.sims_mult.csv',
    shell:
        """
        python workflow/scripts/analysis/topo/run_pair_sim.py \
        -p {input} \
        -t {output.stats} \
        -s {output.sims}
        """