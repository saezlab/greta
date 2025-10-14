thr_r2 = 0.1

rule mdl_o_pearson:
    threads: 16
    singularity: 'workflow/envs/gretabench.sif'
    input:
        mdata=rules.extract_case.output.mdata,
        tf=rules.gen_tfs_lambert.output,
        proms=rules.cre_promoters.output,
    output: out='dts/{dat}/cases/{case}/runs/o_pearson.o_pearson.o_pearson.o_pearson.mdl.csv'
    params:
        mode=config['methods']['pearson']['mode'],
        thr_r2=config['methods']['pearson']['thr_r2'],
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    shell:
        """
        python workflow/scripts/mth/correlation.py \
        -i {input.mdata} \
        -t {input.tf} \
        -p {input.proms} \
        -m {params.mode} \
        -r {params.thr_r2} \
        -o {output.out}
        """

rule mdl_o_spearman:
    threads: 16
    singularity: 'workflow/envs/gretabench.sif'
    input:
        mdata=rules.extract_case.output.mdata,
        tf=rules.gen_tfs_lambert.output,
        proms=rules.cre_promoters.output,
    output: out='dts/{dat}/cases/{case}/runs/o_spearman.o_spearman.o_spearman.o_spearman.mdl.csv'
    params:
        mode=config['methods']['spearman']['mode'],
        thr_r2=config['methods']['spearman']['thr_r2'],
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    shell:
        """
        python workflow/scripts/mth/correlation.py \
        -i {input.mdata} \
        -t {input.tf} \
        -p {input.proms} \
        -m {params.mode} \
        -r {params.thr_r2} \
        -o {output.out}
        """
