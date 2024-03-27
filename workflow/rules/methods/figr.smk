rule pre_figr:
    input:
        'datasets/{dataset}/cases/{case}/mdata.h5mu'
    singularity:
        'workflow/envs/figr.sif'
    benchmark:
        'benchmarks/{dataset}.{case}.figr.pre.txt'
    output:
        'datasets/{dataset}/cases/{case}/runs/figr.pre.h5mu'
    params:
        k=10,
    resources:
        mem_mb=64000,
    shell:
        """
        cp {input} {output}
        Rscript workflow/scripts/methods/figr/pre.R \
        {output} \
        {params.k}
        """

rule p2g_figr:
    input:
        'datasets/{dataset}/cases/{case}/runs/{pre}.pre.h5mu'
    singularity:
        'workflow/envs/figr.sif'
    benchmark:
        'benchmarks/{dataset}.{case}.{pre}.figr.p2g.txt'
    output:
        'datasets/{dataset}/cases/{case}/runs/{pre}.figr.p2g.csv'
    params:
        organism=lambda w: config['datasets'][w.dataset]['organism'],
        ext=500000,
        ncres=3,  # TODO: change to 10
    shell:
        """
        Rscript workflow/scripts/methods/figr/p2g.R \
        {input} \
        {params.organism} \
        {params.ext} \
        {params.ncres} \
        {output}
        """