rule pre_figr:
    input:
        'datasets/{dataset}/cases/{case}/mdata.h5mu'
    singularity:
        'workflow/envs/figr.sif'
    benchmark:
        'benchmarks/{dataset}.{case}.figr.pre.txt'
    output:
        'datasets/{dataset}/cases/{case}/runs/figr.pre.h5mu'
    resources:
        mem_mb=128000,
    shell:
        """
        cp {input} {output}
        Rscript workflow/scripts/methods/figr/pre.R \
        {output}
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
        thr_p2g_pval=0.1,
        ncres=3,  # TODO: change to 10
    resources:
        mem_mb=128000,
        runtime=360,
    shell:
        """
        Rscript workflow/scripts/methods/figr/p2g.R \
        {input} \
        {params.organism} \
        {params.ext} \
        {params.thr_p2g_pval} \
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
        cellK=10,
        dorcK=3,  # TODO: change to 30
    resources:
        mem_mb=128000,
    shell:
        """
        Rscript workflow/scripts/methods/figr/tfb.R \
        {input.d} \
        {params.organism} \
        {input.p} \
        {params.cellK} \
        {params.dorcK} \
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
        cellK=10,
        thr_score=1,
    resources:
        mem_mb=256000,
    shell:
        """
        Rscript workflow/scripts/methods/figr/mdl.R \
        {input.d} \
        {input.p} \
        {input.t} \
        {params.cellK} \
        {params.thr_score} \
        {output}
        """

rule src_figr:
    input:
        'datasets/{dataset}/cases/{case}/mdata.h5mu',
    singularity:
        'workflow/envs/figr.sif'
    benchmark:
        'benchmarks/{dataset}.{case}.figr.src.txt'
    output:
        'datasets/{dataset}/cases/{case}/runs/figr.src.csv'
    params:
        cellK=10,
        organism=lambda w: config['datasets'][w.dataset]['organism'],
        ext=500000,
        thr_p2g_pval=0.1,
        ncres=3,
        dorcK=3,
        thr_score=1,
    resources:
        mem_mb=128000,
    shell:
        """
        Rscript workflow/scripts/methods/figr/src.R \
        {input} \
        {params.cellK} \
        {params.organism} \
        {params.ext} \
        {params.thr_p2g_pval} \
        {params.ncres} \
        {params.dorcK} \
        {params.thr_score} \
        {output}
        """
