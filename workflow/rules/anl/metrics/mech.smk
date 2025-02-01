localrules: mech_tfa
rule mech_tfa:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        grn=lambda wildcards: rules.grn_run.output.out.format(**wildcards),
        rsc=rules.prt_knocktf.output.dir,
    output:
        out='anl/metrics/mech/tfa/{db}/{dat}.{case}/{pre}.{p2g}.{tfb}.{mdl}.scores.csv'
    shell:
        """
        python workflow/scripts/anl/metrics/mech/tfa.py \
        -i {input.grn} \
        -b {input.rsc} \
        -o {output.out}
        """


rule mech_prt:
    threads: 16
    singularity: 'workflow/envs/gretabench.sif'
    input:
        grn=lambda wildcards: rules.grn_run.output.out.format(**wildcards),
        rsc=rules.prt_knocktf.output.dir,
    output:
        out='anl/metrics/mech/prt/{db}/{dat}.{case}/{pre}.{p2g}.{tfb}.{mdl}.scores.csv'
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'] * 2,
    shell:
        """
        set +e
        timeout $(({resources.runtime}-20))m \
        python workflow/scripts/anl/metrics/mech/prt.py \
        -i {input.grn} \
        -b {input.rsc} \
        -o {output.out}
        if [ $? -eq 124 ]; then
            awk 'BEGIN {{ print "name,prc,rcl,f01" }}' > {output.out}
        fi
        """


rule extract_mech_tfm:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        mdata=rules.extract_case.output.mdata,
        tf=rules.gen_tfs_lambert.output,
    output: 'anl/metrics/mech/sss/sss/{dat}.{case}/tfm.csv'
    shell:
        """
        python workflow/scripts/anl/metrics/mech/tfm.py {input.mdata} {input.tf} {output}
        """


rule mech_sss:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        grn=lambda wildcards: rules.grn_run.output.out.format(**wildcards),
        tfm=rules.extract_mech_tfm.output,
    output:
        out='anl/metrics/mech/sss/sss/{dat}.{case}/{pre}.{p2g}.{tfb}.{mdl}.scores.csv'
    params:
        thr_pval=0.01,
    resources:
        mem_mb=8000,
        runtime=60,
    shell:
        """
        set +e
        timeout $(({resources.runtime}-20))m \
    	python workflow/scripts/anl/metrics/mech/sim.py {input.grn} {input.tfm} {params.thr_pval} {output.out}
        if [ $? -eq 124 ]; then
            awk 'BEGIN {{ print "name,prc,rcl,f01" }}' > {output.out}
        fi
    	"""



