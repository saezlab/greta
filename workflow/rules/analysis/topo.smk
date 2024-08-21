import numpy as np


rule topo_mult:
    threads: 32
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
        expand(['datasets/{{dataset}}/cases/{{case}}/runs/{mth}.grn.csv'], mth=['{m}.{m}.{m}.{m}'.format(m='o_' + m) for m in mthds]) + \
        ['datasets/{dataset}/cases/{case}/runs/random.random.random.random.grn.csv',
        'datasets/{dataset}/cases/{case}/runs/collectri.collectri.collectri.collectri.grn.csv',
        'datasets/{dataset}/cases/{case}/runs/dorothea.dorothea.dorothea.dorothea.grn.csv']
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
        expand(['datasets/{{dataset}}/cases/{{case}}/runs/{mth}.grn.csv'], mth=['{m}.{m}.{m}.{m}'.format(m='o_' + m) for m in mthds]) + \
        expand(['datasets/{{dataset}}/cases/{{case}}/runs/{pre}.{p2g}.{tfb}.{mdl}.grn.csv'], zip, pre=mthds, p2g=mthds, tfb=mthds, mdl=mthds) + \
        ['datasets/{dataset}/cases/{case}/runs/random.random.random.random.grn.csv',
        'datasets/{dataset}/cases/{case}/runs/collectri.collectri.collectri.collectri.grn.csv',
        'datasets/{dataset}/cases/{case}/runs/dorothea.dorothea.dorothea.dorothea.grn.csv']
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
