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
