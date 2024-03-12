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
        'benchmarks/{dataset}.{case}.pando.{p2g}.{tfb}.{mdl}.pre.txt'
    output:
        p='datasets/{dataset}/cases/{case}/runs/pando.{p2g}.{tfb}.{mdl}.peaks.csv',
        d='datasets/{dataset}/cases/{case}/runs/pando.{p2g}.{tfb}.{mdl}.pre.h5mu'
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
        d='datasets/{dataset}/cases/{case}/runs/{pre}.pando.{tfb}.{mdl}.pre.h5mu',
        h='gdata/granges/hg38_ensdb_v86.csv',
        m='gdata/granges/mm10_ensdb_v79.csv',
    singularity:
        'workflow/envs/pando.sif'
    benchmark:
        'benchmarks/{dataset}.{case}.{pre}.pando.{tfb}.{mdl}.pre.txt'
    output:
        'datasets/{dataset}/cases/{case}/runs/{pre}.pando.{tfb}.{mdl}.p2g.csv'
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
        d='datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.pando.{mdl}.pre.h5mu',
        p='datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.pando.{mdl}.p2g.csv'
    singularity:
        'workflow/envs/pando.sif'
    benchmark:
        'benchmarks/{dataset}.{case}.{pre}.{p2g}.pando.{mdl}.pre.txt'
    output:
        'datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.pando.{mdl}.tfb.csv'
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
