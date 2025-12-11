rule mdl_o_directnet:
    threads: 1
    singularity: 'workflow/envs/directnet.sif'
    input:
        img='workflow/envs/directnet.sif',
        mdata=rules.extract_case.output.mdata,
        ann=rules.gen_ann_pando.output,
        tss=rules.gen_tss_directnet.output,
    output:
        out='dts/{org}/{dat}/cases/{case}/runs/o_directnet.o_directnet.o_directnet.o_directnet.mdl.csv'
    params:
        ext=config['methods']['directnet']['ext'] // 2,
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'] * 2,
    shell:
        """
        Rscript workflow/scripts/mth/directnet/src.R \
        {input.mdata} \
        {input.ann} \
        {input.tss} \
        {params.ext} \
        {output.out}
        """
