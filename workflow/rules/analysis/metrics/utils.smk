rule aggr_type_task_resource:
    input:
        make_combs(
            path='analysis/metrics/{type}/{task}/{resource}/{dataset}.{case}',
            mthds=mthds,
            name='scores',
            add_src=False
        )
    output:
        'analysis/metrics/{type}/{task}/{resource}/{dataset}.{case}.scores.csv'
    shell:
        """
        python workflow/scripts/analysis/metrics/aggregate.py \
        -i {input} \
        -o {output}
        """

f_scores = 

rule summary:
    input:
        [
            'analysis/mech/prtrb/knocktf/{dataset}.{case}.scores.csv',
            'analysis/mech/tfact/knocktf/{dataset}.{case}.scores.csv',
            'analysis/pred/omics/gtf/{dataset}.{case}.scores.csv',
            'analysis/pred/omics/cretf/{dataset}.{case}.scores.csv',
            'analysis/pred/omics/gcre/{dataset}.{case}.scores.csv',
            'analysis/pred/pathway/kegg/{dataset}.{case}.scores.csv',
            'analysis/pred/pathway/hall/{dataset}.{case}.scores.csv',
            'analysis/pred/pathway/reac/{dataset}.{case}.scores.csv',
            'analysis/prior/tfm/hpa/{dataset}.{case}.scores.csv',
            'analysis/prior/tfm/tfmdb/{dataset}.{case}.scores.csv',
            'analysis/prior/tfbind/chipatlas/{dataset}.{case}.scores.csv',
            'analysis/prior/tfbind/remap2022/{dataset}.{case}.scores.csv',
            'analysis/prior/tfbind/unibind/{dataset}.{case}.scores.csv',
            'analysis/prior/cre/encode/{dataset}.{case}.scores.csv',
            'analysis/prior/cre/gwascatalogue/{dataset}.{case}.scores.csv',
            'analysis/prior/cre/phastcons/{dataset}.{case}.scores.csv',
            'analysis/prior/cre/zhang21/{dataset}.{case}.scores.csv',
            'analysis/prior/cre/promoters/{dataset}.{case}.scores.csv',
            'analysis/prior/eqtl/eqtlcatalogue/{dataset}.{case}.scores.csv',
        ]
    output:
        'analysis/metrics/summary/{dataset}.{case}.csv'
    shell:
        """
        python path/to/script \
        -i {input} \
        -o {output}
        """
