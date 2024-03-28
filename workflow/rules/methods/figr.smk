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

rule tfb_figr:
    input:
        d='datasets/{dataset}/cases/{case}/runs/{pre}.pre.h5mu',
        p='datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.p2g.csv',
    singularity:
        'workflow/envs/figr.sif'
    benchmark:
        'benchmarks/{dataset}.{case}.{pre}.{p2g}.figr.tfb.txt'
    output:
        'datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.figr.tfb.csv'
    params:
        organism=lambda w: config['datasets'][w.dataset]['organism'],
        k=3  # TODO: change to 30
    shell:
        """
        Rscript workflow/scripts/methods/figr/tfb.R \
        {input.d} \
        {params.organism} \
        {input.p} \
        {params.k} \
        {output}
        """

rule mdl_figr:
    input:
        d='datasets/{dataset}/cases/{case}/runs/{pre}.pre.h5mu',
        p='datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.p2g.csv',
        t='datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.tfb.csv',
    singularity:
        'workflow/envs/figr.sif'
    benchmark:
        'benchmarks/{dataset}.{case}.{pre}.{p2g}.{tfb}.figr.mdl.txt'
    output:
        'datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.figr.mdl.csv'
    params:
        thr=0.75,
    shell:
        """
        Rscript workflow/scripts/methods/figr/mdl.R \
        {input.d} \
        {input.p} \
        {input.t} \
        {params.thr} \
        {output}
        """
