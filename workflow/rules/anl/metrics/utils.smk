localrules: aggr_metric, metric_summ


rule aggr_metric:
    threads: 1
    input:
        lambda w: make_combs_rules(w=w, rule_name='{typ}_{tsk}'.format(typ=w.type, tsk=w.task))
    output:
        'anl/metrics/{type}/{task}/{db}/{org}.{dat}.{case}.scores.csv'
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
            'anl/metrics/genom/tfp/europmc/{dat}.{case}.scores.csv',
            'anl/metrics/genom/tfp/intact/{dat}.{case}.scores.csv',
            'anl/metrics/genom/tfb/chipatlas/{dat}.{case}.scores.csv',
            'anl/metrics/genom/tfb/remap2022/{dat}.{case}.scores.csv',
            'anl/metrics/genom/tfb/unibind/{dat}.{case}.scores.csv',
            'anl/metrics/genom/cre/blacklist/{dat}.{case}.scores.csv',
            'anl/metrics/genom/cre/encode/{dat}.{case}.scores.csv',
            'anl/metrics/genom/cre/gwascatalogue/{dat}.{case}.scores.csv',
            'anl/metrics/genom/cre/phastcons/{dat}.{case}.scores.csv',
            'anl/metrics/genom/cre/zhang21/{dat}.{case}.scores.csv',
            'anl/metrics/genom/cre/promoters/{dat}.{case}.scores.csv',
            'anl/metrics/genom/c2g/eqtlcatalogue/{dat}.{case}.scores.csv',
        ]
    output: 'anl/metrics/summary/{dat}.{case}.csv'
    shell:
        """
        python workflow/scripts/anl/metrics/test.py -m {input} -o {output}
        """
