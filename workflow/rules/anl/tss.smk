localrules: tss_aggr


rule tss_gocoef:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        tss_a='dbs/hg38/gen/tss/{mth_a}.bed',
        tss_b='dbs/hg38/gen/tss/{mth_b}.bed',
    output: temp(local('anl/tss/ocoef/{mth_a}.{mth_b}.csv'))
    resources:
        mem_mb=2000,
        runtime=config['max_mins_per_step'],
    shell:
        """
        python workflow/scripts/anl/tss/gocoef.py \
        -a {input.tss_a} \
        -b {input.tss_b} \
        -o {output}
        """


tss_paths = [f'anl/tss/ocoef/{mth_a}.{mth_b}.csv' for mth_a, mth_b in combinations([x for x in mthds + baselines], 2)]
rule tss_aggr:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: tss_paths
    output: "anl/tss/ocoef.csv"
    shell:
        """
        python -c "import pandas as pd; import sys; \
        tss_paths = sys.argv[1:]; \
        df = pd.concat([pd.read_csv(tss_path) for tss_path in tss_paths]); \
        df.to_csv('{output}', index=False);" {input}
        """


rule tss_dist:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        c=rules.tss_aggr.output,
        g='anl/topo/{dat}.{case}.stats_mult.csv'
    output: "anl/tss/{dat}.{case}.dist.csv"
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    params:
        b=baselines,
    shell:
        """
        python workflow/scripts/anl/tss/dist.py \
        -g {input.g} \
        -b {params.b} \
        -o {output}
        """
