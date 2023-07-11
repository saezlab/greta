rule run_pando:
    input:
        data="resources/{dataset}/{trajectory}/mdata.h5mu"
    singularity:
        "workflow/envs/pando.sif"
    benchmark:
        "benchmarks/pando/{dataset}.{trajectory}.run_pando.txt"
    output:
        grn="resources/{dataset}/{trajectory}/pando/grn.csv",
        tri="resources/{dataset}/{trajectory}/pando/tri.csv"
    params:
        organism=lambda w: config[w.dataset]['organism'],
        p_thresh=lambda w: config[w.dataset]['trajectories'][w.trajectory]['pando']['p_thresh'],
        rsq_thresh=lambda w: config[w.dataset]['trajectories'][w.trajectory]['pando']['rsq_thresh'],
        nvar_thresh=lambda w: config[w.dataset]['trajectories'][w.trajectory]['pando']['nvar_thresh']
    shell:
        """
        Rscript workflow/scripts/pando/run_pando.R {input.data} {params.organism} {params.p_thresh} {params.rsq_thresh} {params.nvar_thresh} {output.grn} {output.tri}
        """
