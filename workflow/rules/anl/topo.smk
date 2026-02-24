localrules: topo_diff_seed, topo_inter, topo_fvsd, topo_simul, topo_aggr

n_seeds = 10
rule topo_diff_seed:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: [['dts/{org}/{dts}/cases/16384_16384_' + f'{i}/runs/o_{mth}.o_{mth}.o_{mth}.o_{mth}.grn.csv' for mth in mthds] for i in range(n_seeds)]
    output: 'anl/topo/{org}.{dts}.topo_diff_seed.csv'
    shell:
        """
        python workflow/scripts/anl/topo/stab_seed.py {output}
        """

rule topo_mult:
    threads: 4
    singularity: 'workflow/envs/gretabench.sif'
    input:
        lambda w: make_combs_rules(w=w, rule_name='grn_run')
    output:
        stats='anl/topo/{org}.{dat}.{case}.stats_mult.csv',
        sims='anl/topo/{org}.{dat}.{case}.sims_mult.csv',
    resources:
        mem_mb=128000
    shell:
        """
        python workflow/scripts/anl/topo/run_pair_sim.py \
        -t {output.stats} \
        -s {output.sims}
        """

def find_topo_aggr_paths(dat):
    org = config['dts'][dat]['organism']
    case = 'all'
    return f'anl/topo/{org}.{dat}.{case}.stats_mult.csv'

rule topo_aggr:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: [find_topo_aggr_paths(dat=d) for d in config['dts']]
    output: 'anl/topo/stats.csv'
    shell:
        """
        python workflow/scripts/anl/topo/aggr_stats.py \
        -i {input} \
        -o {output}
        """


rule topo_fvsd:
    threads: 4
    singularity: 'workflow/envs/gretabench.sif'
    input:
        stats=rules.topo_mult.output.stats,
        sims=rules.topo_mult.output.sims,
    output: 'anl/topo/{org}.{dat}.{case}.fvsd.csv',
    shell:
        """
        python workflow/scripts/anl/topo/fvsd.py {input.sims} {input.stats} {output}
        """


rule topo_inter:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: lambda w: make_combs_rules(w=w, rule_name='grn_run')
    output: 'anl/topo/{org}.{dat}.{case}.inter.csv',
    params: min_prop=config['topo_min_prop']
    shell:
        """
        python workflow/scripts/anl/topo/inter.py \
        -g {input} \
        -b {baselines} \
        -p {params.min_prop} \
        -o {output}
        """


rule topo_simul:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: expand('dts/sim/seed_1/{sim_mth}.csv', sim_mth=['celloracle', 'figr', 'pando', 'pearson', 'spearman', 'grnboost', 'random'])
    output: 'anl/topo/sim.1.ocoeff.csv'
    shell:
        """
        python workflow/scripts/anl/topo/simul.py {output}
        """

