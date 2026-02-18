localrules: fig_eval, fig_dts_qc


# Datasets to include in dts_qc figure (excluding fakepitupair and pitunpair)
DTS_QC_DATASETS = [
    ('hg38', 'brain'),
    ('hg38', 'breast'),
    ('hg38', 'embryo'),
    ('hg38', 'eye'),
    ('hg38', 'heart'),
    ('hg38', 'kidney'),
    ('hg38', 'lung'),
    ('hg38', 'pbmc10k'),
    ('hg38', 'pitupair'),
    ('hg38', 'reprofibro'),
    ('hg38', 'skin'),
    ('mm10', 'epalate'),
]


rule fig_dts_qc:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        mdata=expand('dts/{org}/{dat}/cases/all/mdata.h5mu',
            zip, org=[d[0] for d in DTS_QC_DATASETS], dat=[d[1] for d in DTS_QC_DATASETS]),
        qc=expand('anl/dts/{org}.{dat}.all.qc.csv',
            zip, org=[d[0] for d in DTS_QC_DATASETS], dat=[d[1] for d in DTS_QC_DATASETS]),
        nc=expand('anl/dts/{org}.{dat}.all.nc.csv',
            zip, org=[d[0] for d in DTS_QC_DATASETS], dat=[d[1] for d in DTS_QC_DATASETS]),
    output: 'plt/eval/dts_qc.png'
    shell:
        """
        python workflow/scripts/plt/eval/dts_qc.py \
        -c config/config.yaml \
        -o {output}
        """


rule fig_eval:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        metrics=rules.metric_summ.output.metrics,
        scale=rules.metric_summ.output.scale,
        pair=rules.metric_summ.output.pair,
        stats=rules.topo_aggr.output,
        ndbs=rules.dbs_n_per_dts.output,
        stab_seed='anl/stab/unsmthds/hg38.pitupair.scores.csv',
    output: 'plt/eval/fig.pdf'
    shell:
        """
        python workflow/scripts/plt/eval/ranking_figure.py \
        -i {input.metrics} \
        -s {input.scale} \
        -p {input.pair} \
        -t {input.stats} \
        -d {input.ndbs} \
        -b {input.stab_seed} \
        -o {output}
        """
