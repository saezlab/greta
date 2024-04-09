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
        d='datasets/{dataset}/cases/{case}/runs/pando.pre.h5mu'
    params:
        organism=lambda w: config['datasets'][w.dataset]['organism'],
        exclude_exons='TRUE',
    shell:
        """
        Rscript workflow/scripts/methods/pando/pre.R \
        {input.d} \
        {params.organism} \
        {input.h} \
        {input.m} \
        {params.exclude_exons} \
        {output.p}
        python workflow/scripts/methods/pando/pre.py \
        -i {input.d} \
        -p {output.p} \
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
        ext=1e6,
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
        mem_mb=64000,
        runtime=180,
    shell:
        """
        Rscript workflow/scripts/methods/pando/mdl.R \
        {input.d} \
        {input.p} \
        {input.t} \
        {output}
        """
