rule topo_mult:
    threads: 32
    input:
        make_combs(
            path='datasets/{dataset}/cases/{case}/runs/',
            mthds=mthds,
            name='grn',
        )
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
