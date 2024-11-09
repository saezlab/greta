localrules: tss_figr, tss_dictys, tss_celloracle, tss_pando, tss_granie, tss_dist, \
tss_collectri, tss_dorothea, tss_random, tss_scenic

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


rule tss_collectri:
    threads: 1
    input: rules.cre_promoters.output[0]
    output: 'gdata/alltss/collectri.bed'
    shell:
        """
        cp {input} {output}
        """


rule tss_dorothea:
    threads: 1
    input: rules.cre_promoters.output[0]
    output: 'gdata/alltss/dorothea.bed'
    shell:
        """
        cp {input} {output}
        """


rule tss_random:
    threads: 1
    input: rules.cre_promoters.output[0]
    output: 'gdata/alltss/random.bed'
    shell:
        """
        cp {input} {output}
        """


rule tss_scenic:
    threads: 1
    input: rules.cre_promoters.output[0]
    output: 'gdata/alltss/scenic.bed'
    shell:
        """
        cp {input} {output}
        """

rule tss_gocoef:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        tss_a='gdata/alltss/{mth_a}.bed',
        tss_b='gdata/alltss/{mth_b}.bed',
    output: temp(local('analysis/tss/ocoef/{mth_a}.{mth_b}.csv'))
    resources:
        mem_mb=2000,
        runtime=config['max_mins_per_step'],
    shell:
        """
        python workflow/scripts/analysis/tss/gocoef.py \
        -a {input.tss_a} \
        -b {input.tss_b} \
        -o {output}
        """


tss_paths = [f'analysis/tss/ocoef/{mth_a}.{mth_b}.csv' for mth_a, mth_b in combinations([x for x in mthds + baselines if x != 'scenicplus'], 2)]
rule tss_aggr:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: tss_paths
    output: "analysis/tss/ocoef.csv"
    resources:
        mem_mb=2000,
        runtime=config['max_mins_per_step'],
    shell:
        """
        python workflow/scripts/analysis/tss/gocoef.py \
        -t {input} \
        -o {output}
        """

rule tss_dist:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        c='analysis/tss/ocoef.csv',
        g='analysis/topo/{dataset}.{case}.stats_mult.csv'
    output: "analysis/tss/{dataset}.{case}.dist.csv"
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    params:
        b=baselines,
    shell:
        """
        python workflow/scripts/analysis/tss/dist.py \
        -g {input.g} \
        -b {params.b} \
        -o {output}
        """
