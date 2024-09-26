localrules: download_lambert


rule download_lambert:
    output:
        out='gdata/tfs/lambert.csv'
    params:
        url='http://humantfs.ccbr.utoronto.ca/download/v_1.01/TF_names_v_1.01.txt'
    resources:
        mem_mb=2000
    shell:
        "wget '{params.url}' -O {output.out}"


rule mdl_random:
    input:
        mdata=rules.extract_case.output.mdata,
        tf=rules.download_lambert.output.out,
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        out='datasets/{dataset}/cases/{case}/runs/random.random.random.random.mdl.csv'
    params:
        g_perc=0.25,
        scale=1,
        tf_g_ratio=0.10,
        seed=42,
    shell:
        """
        python workflow/scripts/methods/random/grn.py \
        -i {input.mdata} \
        -t {input.tf} \
        -g {params.g_perc} \
        -n {params.scale} \
        -r {params.tf_g_ratio} \
        -s {params.seed} \
        -o {output.out}
        """
