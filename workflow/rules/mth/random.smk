rule mdl_o_random:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        mdata=rules.extract_case.output.mdata,
        tf=lambda w: rules.gen_tfs_lambert_mm10.output if config['dts'][w.dat]['organism'] == 'mm10' else rules.gen_tfs_lambert.output,
        cg=lambda w: rules.cre_promoters_mm10.output if config['dts'][w.dat]['organism'] == 'mm10' else rules.cre_promoters.output,
    output: out='dts/{org}/{dat}/cases/{case}/runs/o_random.o_random.o_random.o_random.mdl.csv'
    params:
        g_perc=config['methods']['random']['g_perc'],
        scale=config['methods']['random']['scale'],
        tf_g_ratio=config['methods']['random']['tf_g_ratio'],
        w_size=config['methods']['random']['w_size'] // 2,
        seed=config['methods']['random']['seed'],
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
