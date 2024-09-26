localrules: aggr_type_task_resource


rule aggr_type_task_resource:
    input:
        lambda w: make_combs_rules(w=w, mthds=mthds, rule_name='{typ}_{tsk}'.format(typ=w.type, tsk=w.task))
    output:
        'analysis/metrics/{type}/{task}/{resource}/{dataset}.{case}.scores.csv'
    shell:
        """
        python workflow/scripts/analysis/metrics/aggregate.py \
        -i {input} \
        -o {output}
        """


rule summary:
    input:
        [
            'analysis/metrics/mech/prtrb/knocktf/{dataset}.{case}.scores.csv',
            'analysis/metrics/mech/tfact/knocktf/{dataset}.{case}.scores.csv',
            'analysis/metrics/pred/omics/gtf/{dataset}.{case}.scores.csv',
            'analysis/metrics/pred/omics/cretf/{dataset}.{case}.scores.csv',
            'analysis/metrics/pred/omics/gcre/{dataset}.{case}.scores.csv',
            'analysis/metrics/pred/pathway/kegg/{dataset}.{case}.scores.csv',
            'analysis/metrics/pred/pathway/hall/{dataset}.{case}.scores.csv',
            'analysis/metrics/pred/pathway/reac/{dataset}.{case}.scores.csv',
            'analysis/metrics/pred/pathway/progeny/{dataset}.{case}.scores.csv',
            'analysis/metrics/prior/tfm/hpa/{dataset}.{case}.scores.csv',
            'analysis/metrics/prior/tfm/tfmdb/{dataset}.{case}.scores.csv',
            'analysis/metrics/prior/tfbind/chipatlas/{dataset}.{case}.scores.csv',
            'analysis/metrics/prior/tfbind/remap2022/{dataset}.{case}.scores.csv',
            'analysis/metrics/prior/tfbind/unibind/{dataset}.{case}.scores.csv',
            'analysis/metrics/prior/cre/encode/{dataset}.{case}.scores.csv',
            'analysis/metrics/prior/cre/gwascatalogue/{dataset}.{case}.scores.csv',
            'analysis/metrics/prior/cre/phastcons/{dataset}.{case}.scores.csv',
            'analysis/metrics/prior/cre/zhang21/{dataset}.{case}.scores.csv',
            'analysis/metrics/prior/cre/promoters/{dataset}.{case}.scores.csv',
            'analysis/metrics/prior/eqtl/eqtlcatalogue/{dataset}.{case}.scores.csv',
        ]
    output:
        'analysis/metrics/summary/{dataset}.{case}.csv'
    shell:
        """
        echo 'Done'
        """
