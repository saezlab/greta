rule run_pando:
    input:
        data="resources/{dataset}/{trajectory}/mdata.h5mu"
    singularity:
        "envs/pando.sif"
    benchmark:
        "benchmarks/pando/{dataset}.{trajectory}.run_pando.txt"
    output:
        grn="resources/{dataset}/{trajectory}/pando/grn.csv",
        tri="resources/{dataset}/{trajectory}/pando/tri.csv"
    params:
        organism=lambda w: config[w.dataset]['organism'],
        p_thresh=lambda w: config[w.dataset]['trajectories'][w.trajectory]['pando']['p_thresh'],
        rsq_thresh=lambda w: config[w.dataset]['trajectories'][w.trajectory]['pando']['rsq_thresh'],
        nvar_thresh=lambda w: config[w.dataset]['trajectories'][w.trajectory]['pando']['nvar_thresh'],
        exclude_exons=lambda w: config[w.dataset]['trajectories'][w.trajectory]['pando']['exclude_exons']
    shell:
        """
        Rscript scripts/pando/run_pando.R {input.data} {params.organism} {params.p_thresh} {params.rsq_thresh} {params.nvar_thresh} {params.exclude_exons} {output.grn} {output.tri}
        """
