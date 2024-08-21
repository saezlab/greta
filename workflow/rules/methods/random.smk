rule download_lambert:
    output:
        'gdata/tfs/lambert.csv'
    params:
        url='http://humantfs.ccbr.utoronto.ca/download/v_1.01/TF_names_v_1.01.txt'
    resources:
        mem_mb=2000
    shell:
        "wget '{params.url}' -O {output}"


rule grn_random:
    input:
        data='datasets/{dataset}/cases/{case}/mdata.h5mu',
        tf='gdata/tfs/lambert.csv',
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        'datasets/{dataset}/cases/{case}/runs/random.random.random.random.grn.csv'
    params:
        g_perc=0.25,
        scale=1,
        tf_g_ratio=0.10,
        seed=42,
    shell:
        """
        python workflow/scripts/methods/random/grn.py \
        -i {input.data} \
        -t {input.tf} \
        -g {params.g_perc} \
        -n {params.scale} \
        -r {params.tf_g_ratio} \
        -s {params.seed} \
        -o {output}
        """
