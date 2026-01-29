localrules: aggr_metric, metric_summ


rule aggr_metric:
    threads: 1
    input:
        lambda w: make_combs_rules(w=w, rule_name='{typ}_{tsk}'.format(typ=w.type, tsk=w.task), do_decoupling=False)
    output:
        'anl/metrics/{type}/{task}/{db}/{org}.{dat}.{case}.scores.csv'
    shell:
        """
        python workflow/scripts/anl/metrics/aggregate.py \
        -i {input} \
        -o {output}
        """


def make_metric_rules(dat):
    org = config['dts'][dat]['organism']
    case = 'all'
    if org == 'hg38':
        return [
            f'anl/metrics/mech/prt/knocktf/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/mech/tfa/knocktf/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/mech/sss/sss/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/pred/omics/gtf/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/pred/omics/cretf/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/pred/omics/gcre/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/pred/gsets/kegg/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/pred/gsets/hall/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/pred/gsets/reac/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/pred/gsets/prog/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/prior/tfm/hpa/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/prior/tfm/tfmdb/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/prior/tfp/europmc/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/prior/tfp/intact/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/prior/grn/collectri/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/genom/tfb/chipatlas/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/genom/tfb/remap2022/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/genom/tfb/unibind/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/genom/cre/blacklist/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/genom/cre/encode/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/genom/cre/gwascatalogue/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/genom/cre/phastcons/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/genom/cre/zhang21/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/genom/cre/promoters/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/genom/c2g/eqtlcatalogue/{org}.{dat}.{case}.scores.csv',
        ]
    elif org == 'mm10':
        return [
            f'anl/metrics/mech/prt/knocktf/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/mech/tfa/knocktf/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/mech/sss/sss/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/pred/omics/gtf/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/pred/omics/cretf/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/pred/omics/gcre/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/pred/gsets/hall/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/pred/gsets/reac/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/pred/gsets/prog/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/prior/grn/collectri/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/genom/tfb/chipatlas/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/genom/tfb/remap2022/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/genom/tfb/unibind/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/genom/cre/blacklist/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/genom/cre/encode/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/genom/cre/phastcons/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/genom/cre/promoters/{org}.{dat}.{case}.scores.csv',
        ]


rule metric_summ:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        aggr=[make_metric_rules(dat=dat) for dat in config['dts'].keys()],
        scale='anl/stab/pitupair.ovc.csv'
    output:
        metrics='anl/metrics/summary/metrics.csv',
        scale='anl/metrics/summary/scalability.csv'
    shell:
        """
        python workflow/scripts/anl/metrics/aggr_all.py {output.metrics} && \
        python workflow/scripts/anl/metrics/scalability.py {input.scale} {output.scale}
        """
