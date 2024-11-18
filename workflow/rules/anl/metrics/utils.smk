localrules: aggr_metric


rule aggr_metric:
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


rule summary:
    input:
        [
            'anl/metrics/mech/prt/knocktf/{dataset}.{case}.scores.csv',
            'anl/metrics/mech/tfa/knocktf/{dataset}.{case}.scores.csv',
            'anl/metrics/pred/omics/gtf/{dataset}.{case}.scores.csv',
            'anl/metrics/pred/omics/cretf/{dataset}.{case}.scores.csv',
            'anl/metrics/pred/omics/gcre/{dataset}.{case}.scores.csv',
            'anl/metrics/pred/gsets/kegg/{dataset}.{case}.scores.csv',
            'anl/metrics/pred/gsets/hall/{dataset}.{case}.scores.csv',
            'anl/metrics/pred/gsets/reac/{dataset}.{case}.scores.csv',
            'anl/metrics/pred/gsets/progeny/{dataset}.{case}.scores.csv',
            'anl/metrics/prior/tfm/hpa/{dataset}.{case}.scores.csv',
            'anl/metrics/prior/tfm/tfmdb/{dataset}.{case}.scores.csv',
            'anl/metrics/prior/tfb/chipatlas/{dataset}.{case}.scores.csv',
            'anl/metrics/prior/tfb/remap2022/{dataset}.{case}.scores.csv',
            'anl/metrics/prior/tfb/unibind/{dataset}.{case}.scores.csv',
            'anl/metrics/prior/cre/encode/{dataset}.{case}.scores.csv',
            'anl/metrics/prior/cre/gwascatalogue/{dataset}.{case}.scores.csv',
            'anl/metrics/prior/cre/phastcons/{dataset}.{case}.scores.csv',
            'anl/metrics/prior/cre/zhang21/{dataset}.{case}.scores.csv',
            'anl/metrics/prior/cre/promoters/{dataset}.{case}.scores.csv',
            'anl/metrics/prior/c2g/eqtlcatalogue/{dataset}.{case}.scores.csv',
        ]
    output:
        'anl/metrics/summary/{dat}.{case}.csv'
    shell:
        """
        echo 'Done'
        """
