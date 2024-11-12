localrules: gen_tss_celloracle, gen_tss_dictys, gen_tss_figr, gen_tss_granie, gen_tss_pando, \
gen_tss_collectri, gen_tss_dorothea, gen_tss_random, gen_tss_scenic


rule gen_tss_celloracle:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    output: 'dbs/hg38/gen/tss/celloracle.bed'
    shell:
        """
        python workflow/scripts/dbs/gen/tss/celloracle.py -o {output}
        """


rule gen_tss_dictys:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: rules.gen_ann_dictys.output
    output: 'dbs/hg38/gen/tss/dictys.bed'
    shell:
        """
        python workflow/scripts/dbs/gen/tss/dictys.py \
        -i {input} \
        -o {output}
        """


rule gen_tss_figr:
    threads: 1
    singularity: 'workflow/envs/figr.sif'
    output: 'dbs/hg38/gen/tss/figr.bed'
    shell:
        """
        Rscript workflow/scripts/dbs/gen/tss/figr.R {output}
        """


rule gen_tss_pando:
    threads: 1
    singularity: 'workflow/envs/pando.sif'
    output: 'dbs/hg38/gen/tss/pando.bed'
    shell:
        """
        Rscript workflow/scripts/dbs/gen/tss/pando.R {output}
        """


rule gen_tss_granie:
    threads: 1
    singularity: 'workflow/envs/granie.sif'
    output: 'dbs/hg38/gen/tss/granie.bed'
    shell:
        """
        Rscript workflow/scripts/dbs/gen/tss/granie.R {output}
        """


rule gen_tss_collectri:
    threads: 1
    input: rules.cre_promoters.output
    output: 'dbs/hg38/gen/tss/collectri.bed'
    shell:
        """
        cp {input} {output}
        """


rule gen_tss_dorothea:
    threads: 1
    input: rules.cre_promoters.output
    output: 'dbs/hg38/gen/tss/dorothea.bed'
    shell:
        """
        cp {input} {output}
        """


rule gen_tss_random:
    threads: 1
    input: rules.cre_promoters.output
    output: 'dbs/hg38/gen/tss/random.bed'
    shell:
        """
        cp {input} {output}
        """


rule gen_tss_scenic:
    threads: 1
    input: rules.cre_promoters.output
    output: 'dbs/hg38/gen/tss/scenic.bed'
    shell:
        """
        cp {input} {output}
        """
