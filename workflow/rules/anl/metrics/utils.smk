localrules: aggr_type_task_resource


rule aggr_type_task_resource:
    input:
        lambda w: make_combs_rules(w=w, mthds=mthds, rule_name='{typ}_{tsk}'.format(typ=w.type, tsk=w.task))
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
            'anl/metrics/mech/prtrb/knocktf/{dataset}.{case}.scores.csv',
            'anl/metrics/mech/tfact/knocktf/{dataset}.{case}.scores.csv',
            'anl/metrics/pred/omics/gtf/{dataset}.{case}.scores.csv',
            'anl/metrics/pred/omics/cretf/{dataset}.{case}.scores.csv',
            'anl/metrics/pred/omics/gcre/{dataset}.{case}.scores.csv',
            'anl/metrics/pred/pathway/kegg/{dataset}.{case}.scores.csv',
            'anl/metrics/pred/pathway/hall/{dataset}.{case}.scores.csv',
            'anl/metrics/pred/pathway/reac/{dataset}.{case}.scores.csv',
            'anl/metrics/pred/pathway/progeny/{dataset}.{case}.scores.csv',
            'anl/metrics/prior/tfm/hpa/{dataset}.{case}.scores.csv',
            'anl/metrics/prior/tfm/tfmdb/{dataset}.{case}.scores.csv',
            'anl/metrics/prior/tfbind/chipatlas/{dataset}.{case}.scores.csv',
            'anl/metrics/prior/tfbind/remap2022/{dataset}.{case}.scores.csv',
            'anl/metrics/prior/tfbind/unibind/{dataset}.{case}.scores.csv',
            'anl/metrics/prior/cre/encode/{dataset}.{case}.scores.csv',
            'anl/metrics/prior/cre/gwascatalogue/{dataset}.{case}.scores.csv',
            'anl/metrics/prior/cre/phastcons/{dataset}.{case}.scores.csv',
            'anl/metrics/prior/cre/zhang21/{dataset}.{case}.scores.csv',
            'anl/metrics/prior/cre/promoters/{dataset}.{case}.scores.csv',
            'anl/metrics/prior/eqtl/eqtlcatalogue/{dataset}.{case}.scores.csv',
        ]
    output:
        'anl/metrics/summary/{dat}.{case}.csv'
    shell:
        """
        echo 'Done'
        """
