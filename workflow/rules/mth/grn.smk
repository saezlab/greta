localrules: grn_run


rule grn_run:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: lambda wildcards: map_rules('mdl', wildcards.mdl),
    output:
        out='dts/{org}/{dat}/cases/{case}/runs/{pre}.{p2g}.{tfb}.{mdl}.grn.csv'
    shell:
        """
        python workflow/scripts/mth/grn.py \
        -i {input} \
        -o {output.out}
        """


rule mdl_o_collectri:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        mdata=rules.extract_case.output.mdata,
        grn=lambda w: rules.gst_collectri_mm10.output if config['dts'][w.dat]['organism'] == 'mm10' else rules.gst_collectri.output,
        proms=lambda w: rules.cre_promoters_mm10.output if config['dts'][w.dat]['organism'] == 'mm10' else rules.cre_promoters.output,
    output:
        out='dts/{org}/{dat}/cases/{case}/runs/o_collectri.o_collectri.o_collectri.o_collectri.mdl.csv'
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    shell:
        """
        python workflow/scripts/mth/prc_prior_grn.py \
        -g {input.grn} \
        -d {input.mdata} \
        -p {input.proms} \
        -o {output.out}
        """


rule mdl_o_dorothea:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        mdata=rules.extract_case.output.mdata,
        grn=lambda w: rules.gst_dorothea_mm10.output if config['dts'][w.dat]['organism'] == 'mm10' else rules.gst_dorothea.output,
        proms=lambda w: rules.cre_promoters_mm10.output if config['dts'][w.dat]['organism'] == 'mm10' else rules.cre_promoters.output,
    output:
        out='dts/{org}/{dat}/cases/{case}/runs/o_dorothea.o_dorothea.o_dorothea.o_dorothea.mdl.csv'
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    shell:
        """
        python workflow/scripts/mth/prc_prior_grn.py \
        -g {input.grn} \
        -d {input.mdata} \
        -p {input.proms} \
        -o {output.out}
        """
