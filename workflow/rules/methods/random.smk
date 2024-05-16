rule pre_random:
    input:
        'datasets/{dataset}/cases/{case}/mdata.h5mu'
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        'datasets/{dataset}/cases/{case}/runs/random.pre.h5mu'
    shell:
        """
        python workflow/scripts/methods/random/pre.py \
        -i {input} \
        -o {output}
        """

rule p2g_random:
    input:
        'datasets/{dataset}/cases/{case}/runs/{pre}.pre.h5mu'
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        'datasets/{dataset}/cases/{case}/runs/{pre}.random.p2g.csv'
    shell:
        """
        python workflow/scripts/methods/random/p2g.py -i {input} -o {output}
        """

rule download_lambert:
    output:
        'gdata/tfs/lambert.csv'
    params:
        url='http://humantfs.ccbr.utoronto.ca/download/v_1.01/TF_names_v_1.01.txt'
    resources:
        mem_mb=2000
    shell:
        "wget '{params.url}' -O {output}"

rule tfb_random:
    input:
        data='datasets/{dataset}/cases/{case}/runs/{pre}.pre.h5mu',
        p2g='datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.p2g.csv',
        tfs='gdata/tfs/lambert.csv',
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        'datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.random.tfb.csv'
    shell:
        """
        python workflow/scripts/methods/random/tfb.py \
        -i {input.data} \
        -t {input.tfs} \
        -p {input.p2g} \
        -o {output}
        """

rule mdl_random:
    input:
        p2g='datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.p2g.csv',
        tfb='datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.tfb.csv',
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        'datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.random.mdl.csv'
    shell:
        """
        python workflow/scripts/methods/random/mdl.py \
        -p {input.p2g} \
        -t {input.tfb} \
        -o {output}
        """
