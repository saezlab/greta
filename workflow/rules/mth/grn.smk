localrules: grn_run


rule grn_run:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: lambda wildcards: map_rules('mdl', wildcards.mdl),
    output:
        out='dts/{dat}/cases/{case}/runs/{pre}.{p2g}.{tfb}.{mdl}.grn.csv'
    shell:
        """
        python workflow/scripts/mth/grn.py \
        -i {input} \
        -o {output.out}
        """


rule mdl_collectri:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        mdata=rules.extract_case.output.mdata,
        grn=rules.gst_collectri.output,
        proms=rules.cre_promoters.output,
    output:
        out='dts/{dat}/cases/{case}/runs/collectri.collectri.collectri.collectri.mdl.csv'
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


rule mdl_dorothea:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        mdata=rules.extract_case.output.mdata,
        grn=rules.gst_dorothea.output,
        proms=rules.cre_promoters.output,
    output:
        out='dts/{dat}/cases/{case}/runs/dorothea.dorothea.dorothea.dorothea.mdl.csv'
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
