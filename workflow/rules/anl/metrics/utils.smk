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
            'anl/metrics/mech/prt/knocktf/{org}.{dat}.{case}.scores.csv',
            'anl/metrics/mech/tfa/knocktf/{org}.{dat}.{case}.scores.csv',
            'anl/metrics/mech/sss/sss/{org}.{dat}.{case}.scores.csv',
            'anl/metrics/pred/omics/gtf/{org}.{dat}.{case}.scores.csv',
            'anl/metrics/pred/omics/cretf/{org}.{dat}.{case}.scores.csv',
            'anl/metrics/pred/omics/gcre/{org}.{dat}.{case}.scores.csv',
            'anl/metrics/pred/gsets/kegg/{org}.{dat}.{case}.scores.csv',
            'anl/metrics/pred/gsets/hall/{org}.{dat}.{case}.scores.csv',
            'anl/metrics/pred/gsets/reac/{org}.{dat}.{case}.scores.csv',
            'anl/metrics/pred/gsets/prog/{org}.{dat}.{case}.scores.csv',
            'anl/metrics/prior/tfm/hpa/{org}.{dat}.{case}.scores.csv',
            'anl/metrics/prior/tfm/tfmdb/{org}.{dat}.{case}.scores.csv',
            'anl/metrics/prior/tfp/europmc/{org}.{dat}.{case}.scores.csv',
            'anl/metrics/prior/tfp/intact/{org}.{dat}.{case}.scores.csv',
            'anl/metrics/prior/grn/collectri/{org}.{dat}.{case}.scores.csv',
            'anl/metrics/genom/tfb/chipatlas/{org}.{dat}.{case}.scores.csv',
            'anl/metrics/genom/tfb/remap2022/{org}.{dat}.{case}.scores.csv',
            'anl/metrics/genom/tfb/unibind/{org}.{dat}.{case}.scores.csv',
            'anl/metrics/genom/cre/blacklist/{org}.{dat}.{case}.scores.csv',
            'anl/metrics/genom/cre/encode/{org}.{dat}.{case}.scores.csv',
            'anl/metrics/genom/cre/gwascatalogue/{org}.{dat}.{case}.scores.csv',
            'anl/metrics/genom/cre/phastcons/{org}.{dat}.{case}.scores.csv',
            'anl/metrics/genom/cre/zhang21/{org}.{dat}.{case}.scores.csv',
            'anl/metrics/genom/cre/promoters/{org}.{dat}.{case}.scores.csv',
            'anl/metrics/genom/c2g/eqtlcatalogue/{org}.{dat}.{case}.scores.csv',
        ]
    output: 'anl/metrics/summary/{org}.{dat}.{case}.csv'
    shell:
        """
        python workflow/scripts/anl/metrics/test.py -m {input} -o {output}
        """
