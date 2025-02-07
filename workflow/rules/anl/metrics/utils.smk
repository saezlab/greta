localrules: aggr_metric, metric_summ


rule aggr_metric:
    threads: 1
    input:
        lambda w: make_combs_rules(w=w, mthds=mthds, baselines=baselines, rule_name='{typ}_{tsk}'.format(typ=w.type, tsk=w.task))
    output:
        'anl/metrics/{type}/{task}/{db}/{dat}.{case}.scores.csv'
    shell:
        """
        python workflow/scripts/anl/metrics/aggregate.py \
        -i {input} \
        -o {output}
        """


rule metric_summ:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        [
            'anl/metrics/mech/prt/knocktf/{dat}.{case}.scores.csv',
            'anl/metrics/mech/tfa/knocktf/{dat}.{case}.scores.csv',
            'anl/metrics/mech/sss/sss/{dat}.{case}.scores.csv',
            'anl/metrics/pred/omics/gtf/{dat}.{case}.scores.csv',
            'anl/metrics/pred/omics/cretf/{dat}.{case}.scores.csv',
            'anl/metrics/pred/omics/gcre/{dat}.{case}.scores.csv',
            'anl/metrics/pred/gsets/kegg/{dat}.{case}.scores.csv',
            'anl/metrics/pred/gsets/hall/{dat}.{case}.scores.csv',
            'anl/metrics/pred/gsets/reac/{dat}.{case}.scores.csv',
            'anl/metrics/pred/gsets/prog/{dat}.{case}.scores.csv',
            'anl/metrics/prior/tfm/hpa/{dat}.{case}.scores.csv',
            'anl/metrics/prior/tfm/tfmdb/{dat}.{case}.scores.csv',
            'anl/metrics/prior/tfp/europmc/{dat}.{case}.scores.csv',
            'anl/metrics/prior/tfp/intact/{dat}.{case}.scores.csv',
            'anl/metrics/prior/tfb/chipatlas/{dat}.{case}.scores.csv',
            'anl/metrics/prior/tfb/remap2022/{dat}.{case}.scores.csv',
            'anl/metrics/prior/tfb/unibind/{dat}.{case}.scores.csv',
            'anl/metrics/prior/cre/blacklist/{dat}.{case}.scores.csv',
            'anl/metrics/prior/cre/encode/{dat}.{case}.scores.csv',
            'anl/metrics/prior/cre/gwascatalogue/{dat}.{case}.scores.csv',
            'anl/metrics/prior/cre/phastcons/{dat}.{case}.scores.csv',
            'anl/metrics/prior/cre/zhang21/{dat}.{case}.scores.csv',
            'anl/metrics/prior/cre/promoters/{dat}.{case}.scores.csv',
            'anl/metrics/prior/c2g/eqtlcatalogue/{dat}.{case}.scores.csv',
        ]
    output: 'anl/metrics/summary/{dat}.{case}.csv'
    shell:
        """
        python workflow/scripts/anl/metrics/test.py -m {input} -o {output}
        """
