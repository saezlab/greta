localrules: net_centr


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

rule net_centr:
    threads: 1
    input:
        grn_paths=lambda w: make_combs_rules(w=w, mthds=mthds, baselines=baselines, rule_name='grn_run')
    output:
        stats='analysis/topo/{dataset}.{case}.net_centr.csv',
    run:
        import pandas as pd
        import igraph as ig
        import numpy as np
        from tqdm import tqdm
        df = []
        names = []
        on = ['source', 'target']
        for grn_path in tqdm(input['grn_paths']):
            name = grn_path.split('.')[-3].replace('o_', '')
            grn = pd.read_csv(grn_path).drop_duplicates(on)
            tfs = set(grn['source']) & set(grn['target'])
            msk = grn['source'].isin(tfs) & grn['target'].isin(tfs)
            grn = grn.loc[msk, :]
            g = ig.Graph.TupleList(list(zip(grn['source'], grn['target'])), directed=True)
            bet = np.mean(g.betweenness())
            odg = np.mean(g.outdegree())
            if not g.is_acyclic():
                eig = np.mean(g.eigenvector_centrality())
            else:
                eig = 0.
            df.append([name, odg, bet, eig])
        df = pd.DataFrame(df, columns=['name', 'odg', 'bet', 'eig'])
        df.to_csv(output['stats'], index=False)
