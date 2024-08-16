rule download_granges:
    output:
        h='gdata/granges/hg38_ensdb_v86.csv',
        m='gdata/granges/mm10_ensdb_v79.csv',
        d=directory('gdata/granges'),
    singularity:
        'workflow/envs/pando.sif'
    shell:
        """
        Rscript workflow/scripts/methods/pando/get_granges.R {output.h} {output.m}
        """

rule pre_pando:
    input:
        d='datasets/{dataset}/cases/{case}/mdata.h5mu',
        h='gdata/granges/hg38_ensdb_v86.csv',
        m='gdata/granges/mm10_ensdb_v79.csv',
    singularity:
        'workflow/envs/pando.sif'
    benchmark:
        'benchmarks/{dataset}.{case}.pando.pre.txt'
    output:
        p=temp(local('datasets/{dataset}/cases/{case}/runs/pando.peaks.csv')),
        m=temp(local('datasets/{dataset}/cases/{case}/runs/pando.matches.csv')),
        d='datasets/{dataset}/cases/{case}/runs/pando.pre.h5mu'
    params:
        organism=lambda w: config['datasets'][w.dataset]['organism'],
        exclude_exons=config['methods']['pando']['exclude_exons'],
    shell:
        """
        Rscript workflow/scripts/methods/pando/pre.R \
        {input.d} \
        {params.organism} \
        {input.h} \
        {input.m} \
        {params.exclude_exons} \
        {output.p} \
        {output.m}
        python workflow/scripts/methods/pando/pre.py \
        -i {input.d} \
        -p {output.p} \
        -m {output.m} \
        -o {output.d}
        """

rule p2g_pando:
    input:
        d='datasets/{dataset}/cases/{case}/runs/{pre}.pre.h5mu',
        h='gdata/granges/hg38_ensdb_v86.csv',
        m='gdata/granges/mm10_ensdb_v79.csv',
    singularity:
        'workflow/envs/pando.sif'
    benchmark:
        'benchmarks/{dataset}.{case}.{pre}.pando.p2g.txt'
    output:
        'datasets/{dataset}/cases/{case}/runs/{pre}.pando.p2g.csv'
    params:
        organism=lambda w: config['datasets'][w.dataset]['organism'],
        ext=config['methods']['pando']['ext'],
    shell:
        """
        Rscript workflow/scripts/methods/pando/p2g.R \
        {input.d} \
        {params.organism} \
        {input.h} \
        {input.m} \
        {params.ext} \
        {output}
        """

rule tfb_pando:
    input:
        d='datasets/{dataset}/cases/{case}/runs/{pre}.pre.h5mu',
        p='datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.p2g.csv'
    singularity:
        'workflow/envs/pando.sif'
    benchmark:
        'benchmarks/{dataset}.{case}.{pre}.{p2g}.pando.tfb.txt'
    output:
        'datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.pando.tfb.csv'
    params:
        organism=lambda w: config['datasets'][w.dataset]['organism']
    shell:
        """
        Rscript workflow/scripts/methods/pando/tfb.R \
        {input.d} \
        {params.organism} \
        {input.p} \
        {output}
        """

rule mdl_pando:
    input:
        d='datasets/{dataset}/cases/{case}/runs/{pre}.pre.h5mu',
        p='datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.p2g.csv',
        t='datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.tfb.csv'
    singularity:
        'workflow/envs/pando.sif'
    benchmark:
        'benchmarks/{dataset}.{case}.{pre}.{p2g}.{tfb}.pando.mdl.txt'
    output:
        'datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.pando.mdl.csv'
    resources:
        mem_mb=512000,
        runtime=1440,
    params:
        thr_corr=config['methods']['pando']['thr_corr'],
        p_thresh=config['methods']['pando']['p_thresh'],
        rsq_thresh=config['methods']['pando']['rsq_thresh'],
        nvar_thresh=config['methods']['pando']['nvar_thresh'],
        min_genes_per_module=config['methods']['pando']['min_genes_per_module'],
    shell:
        """
        Rscript workflow/scripts/methods/pando/mdl.R \
        {input.d} \
        {input.p} \
        {input.t} \
        {params.thr_corr} \
        {params.p_thresh} \
        {params.rsq_thresh} \
        {params.nvar_thresh} \
        {params.min_genes_per_module} \
        {output}
        """

rule src_pando:
    input:
        d='datasets/{dataset}/cases/{case}/mdata.h5mu',
        h='gdata/granges/hg38_ensdb_v86.csv',
        m='gdata/granges/mm10_ensdb_v79.csv',
    singularity:
        'workflow/envs/pando.sif'
    benchmark:
        'benchmarks/{dataset}.{case}.pando.src.txt'
    output:
        'datasets/{dataset}/cases/{case}/runs/pando.src.csv'
    resources:
        mem_mb=512000,
        runtime=720,
    params:
        organism=lambda w: config['datasets'][w.dataset]['organism'],
        exclude_exons=config['methods']['pando']['exclude_exons'],
        ext=config['methods']['pando']['ext'],
        thr_corr=config['methods']['pando']['thr_corr'],
        p_thresh=config['methods']['pando']['p_thresh'],
        rsq_thresh=config['methods']['pando']['rsq_thresh'],
        nvar_thresh=config['methods']['pando']['nvar_thresh'],
        min_genes_per_module=config['methods']['pando']['min_genes_per_module'],
    shell:
        """
        Rscript workflow/scripts/methods/pando/src.R \
        {input.d} \
        {params.exclude_exons} \
        {params.ext} \
        {params.thr_corr} \
        {params.p_thresh} \
        {params.rsq_thresh} \
        {params.nvar_thresh} \
        {params.min_genes_per_module} \
        {params.organism} \
        {input.h} \
        {input.m} \
        {output}
        """
