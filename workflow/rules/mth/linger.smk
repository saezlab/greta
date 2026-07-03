rule mdl_o_linger:
    threads: 24 
    singularity: 'workflow/envs/linger.sif'
    input:
        img='workflow/envs/linger.sif',
        linger_GRN=rules.linger_prior.output.dir,
        mdata=rules.extract_case.output.mdata
    output:
        dir=directory('dts/{org}/{dat}/cases/{case}/runs/linger/'),
        out='dts/{org}/{dat}/cases/{case}/runs/o_linger.o_linger.o_linger.o_linger.mdl.csv'
    params:
        version=config['methods']['linger']['version'],
        mode=config['methods']['linger']['mode'],
        organism=lambda w: config['dts'][w.dat]['organism'],
        script='workflow/scripts/mth/linger/linger.sh'
    resources:
        mem_mb=lambda wildcards, attempt: restart_mem(wildcards, attempt) * 2,
        runtime=(
            30 if config['methods']['linger']['version'] == 'baseline'
            else 140 if config['methods']['linger']['mode'] == 'parallel'
            else 360
        )
    shell:
        """
        mkdir -p {output.dir}
        set -e
        timeout $(({resources.runtime}-20))m \
        bash {params.script} \
        --linger_GRN {input.linger_GRN} \
        --out_dir {output.dir} \
        --path_mdata {input.mdata} \
        --version {params.version} \
        --genome {params.organism} \
        --mode {params.mode} \
        --path_out {output.out}
        if [ $? -eq 124 ]; then
            awk 'BEGIN {{ print "source,target,score,pval" }}' > {output.out}
        fi
        """