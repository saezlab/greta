rule run_pair_sim_mult_runs:
    input:
        expand(['datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.{mdl}.grn.csv'], pre=mthds, p2g=mthds, tfb=mthds, mdl=mthds)
    output:
        'analysis/topo/{dataset}.{case}.mult.csv'
    params:
        m=mthds
    shell:
        """
        python workflow/scripts/analysis/topo/run_pair_sim.py \
        -p {input} \
        -m {params.m} \
        -o {output}
        """

rule run_pair_sim_single_runs:
    input:
        expand(['datasets/{dataset}/cases/{case}/runs/{mth}.src.csv'], mth=mthds) + ['datasets/{dataset}/cases/{case}/runs/random.grn.csv']
    output:
        'analysis/topo/{dataset}.{case}.single.csv'
    params:
        m=mthds
    shell:
        """
        python workflow/scripts/analysis/topo/run_pair_sim.py \
        -p {input} \
        -m {params.m} \
        -o {output}
        """

rule run_pair_sim_orig_runs:
    input:
        expand(['datasets/{dataset}/cases/{case}/runs/{mth}.src.csv'], mth=mthds) + expand(['datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.{mdl}.grn.csv'], zip, pre=mthds, p2g=mthds, tfb=mthds, mdl=mthds)
    output:
        'analysis/topo/{dataset}.{case}.orig.csv'
    params:
        m=mthds
    shell:
        """
        python workflow/scripts/analysis/topo/run_pair_sim.py \
        -p {input} \
        -m {params.m} \
        -o {output}
        """
