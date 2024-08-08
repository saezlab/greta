rule get_grn:
    input:
        'datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.{mdl}.mdl.csv'
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        'datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.{mdl}.grn.csv'
    resources:
        mem_mb=16000,
    shell:
        """
        python workflow/scripts/methods/grn.py \
        -i {input} \
        -o {output}
        """
