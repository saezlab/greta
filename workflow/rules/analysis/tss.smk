localrules: tss_figr, tss_dictys, tss_celloracle, tss_pando, tss_granie

rule tss_figr:
    threads: 1
    singularity: 'workflow/envs/figr.sif'
    output: 'gdata/alltss/figr.bed'
    shell:
        """
        Rscript workflow/scripts/analysis/tss/figr.R {output}
        """


rule tss_dictys:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: rules.download_gene_annotations_dictys.output.bed
    output: 'gdata/alltss/dictys.bed'
    shell:
        """
        python workflow/scripts/analysis/tss/dictys.py \
        -i {input} \
        -o {output}
        """


rule tss_celloracle:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    output: 'gdata/alltss/celloracle.bed'
    shell:
        """
        python workflow/scripts/analysis/tss/celloracle.py -o {output}
        """


rule tss_pando:
    threads: 1
    singularity: 'workflow/envs/pando.sif'
    output: 'gdata/alltss/pando.bed'
    shell:
        """
        Rscript workflow/scripts/analysis/tss/pando.R {output}
        """


rule tss_granie:
    threads: 1
    singularity: 'workflow/envs/granie.sif'
    output: 'gdata/alltss/granie.bed'
    shell:
        """
        Rscript workflow/scripts/analysis/tss/granie.R {output}
        """


rule compare_tss:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: [map_rules('tss', mth, out=0) for mth in mthds if mth != 'scenicplus']  # TODO: remove when ready
    output: "analysis/tss/ocoef.csv"
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    shell:
        """
        python workflow/scripts/analysis/tss/tss.py \
        -t {input} \
        -o {output}
        """
