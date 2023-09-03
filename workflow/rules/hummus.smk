rule run_hummus:
    input:
        data="resources/{dataset}/{trajectory}/mdata.h5mu"
    singularity:
        "envs/hummus.sif"
    benchmark:
        "benchmarks/hummus/{dataset}.{trajectory}.run_hummus.txt"
    output:
        multilayer_f=directory("resources/{dataset}/{trajectory}/hummus/multilayer"),
        grn="resources/{dataset}/{trajectory}/hummus/grn.csv",
        tri="resources/{dataset}/{trajectory}/hummus/tri.csv"
    params:
        genie3_cores = lambda w: config[w.dataset]['trajectories'][w.trajectory]['hummus']['genie3_cores'],
        organism = lambda w: config[w.dataset]['organism'],
        genie3_thresh = lambda w: config[w.dataset]['trajectories'][w.trajectory]['hummus']['genie3_thresh'],
        cicero_k_per_pseudocells = lambda w: config[w.dataset]['trajectories'][w.trajectory]['hummus']['cicero_k_per_pseudocells'],
        cicero_thresh = lambda w: config[w.dataset]['trajectories'][w.trajectory]['hummus']['cicero_thresh'],
        upstream_gene = lambda w: config[w.dataset]['trajectories'][w.trajectory]['hummus']['upstream_gene'],
        downstream_gene = lambda w: config[w.dataset]['trajectories'][w.trajectory]['hummus']['downstream_gene'],
        only_tss = lambda w: config[w.dataset]['trajectories'][w.trajectory]['hummus']['only_tss']
    shell:
        """
        Rscript scripts/hummus/run_hummus.R\
        {input.data}\
        {params.genie3_cores}\
        {params.organism}\
        {params.genie3_thresh}\
        {params.cicero_k_per_pseudocells}\
        {params.cicero_thresh}\
        {params.upstream_gene}\
        {params.downstream_gene}\
        {params.only_tss}\
        {output.multilayer_f}\
        {output.grn}\
        {output.tri}
        """
