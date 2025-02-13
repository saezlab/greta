rule mdl_random:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        mdata=rules.extract_case.output.mdata,
        tf=rules.gen_tfs_lambert.output,
        cg=rules.cre_promoters.output,
    output: out='dts/{dat}/cases/{case}/runs/random.random.random.random.mdl.csv'
    params:
        g_perc=0.25,
        scale=1,
        tf_g_ratio=0.10,
        w_size=250000,
        seed=lambda w: config['dts']['pitupair']['cases'][w.case].get('seed', 42),
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    shell:
        """
        python workflow/scripts/mth/random/grn.py \
        -i {input.mdata} \
        -t {input.tf} \
        -c {input.cg} \
        -g {params.g_perc} \
        -n {params.scale} \
        -r {params.tf_g_ratio} \
        -w {params.w_size} \
        -s {params.seed} \
        -o {output.out}
        """
