import numpy as np


rule topo_mult:
    input:
        expand(['datasets/{{dataset}}/cases/{{case}}/runs/{pre}.{p2g}.{tfb}.{mdl}.grn.csv'], pre=mthds, p2g=mthds, tfb=mthds, mdl=mthds)
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

rule topo_single:
    input:
        expand(['datasets/{{dataset}}/cases/{{case}}/runs/{mth}.src.csv'], mth=mthds) + ['datasets/{dataset}/cases/{case}/runs/random.grn.csv']
    output:
        stats='analysis/topo/{dataset}.{case}.stats_single.csv',
        sims='analysis/topo/{dataset}.{case}.sims_single.csv',
    shell:
        """
        python workflow/scripts/analysis/topo/run_pair_sim.py \
        -p {input} \
        -t {output.stats} \
        -s {output.sims}
        """

rule topo_orign:
    input:
        expand(['datasets/{{dataset}}/cases/{{case}}/runs/{mth}.src.csv'], mth=mthds) + \
        expand(['datasets/{{dataset}}/cases/{{case}}/runs/{pre}.{p2g}.{tfb}.{mdl}.grn.csv'], zip, pre=mthds, p2g=mthds, tfb=mthds, mdl=mthds) + \
        ['datasets/{dataset}/cases/{case}/runs/random.grn.csv']
    output:
        stats='analysis/topo/{dataset}.{case}.stats_orign.csv',
        sims='analysis/topo/{dataset}.{case}.sims_orign.csv',
    shell:
        """
        python workflow/scripts/analysis/topo/run_pair_sim.py \
        -p {input} \
        -t {output.stats} \
        -s {output.sims}
        """
