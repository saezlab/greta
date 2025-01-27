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


rule mech_sim:
    input:
        grn="data/scenic.csv",   #change inputs, add tfs gold dictionary
        data="data/mat.csv"
    output:
        out='anl/metrics/mech/prt/{db}/{dat}.{case}/{pre}.{p2g}.{tfb}.{mdl}.sim.csv'
    shell:
        """
	python workflow/scripts/anl/metrics/mech/sim.py \ 
	-g {input.grn} \
	-d {input.data} \
	-o {output.out} \
	"""



