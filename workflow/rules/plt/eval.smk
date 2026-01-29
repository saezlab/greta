localrules: fig_eval


rule fig_eval:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        metrics=rules.metric_summ.output.metrics,
        scale=rules.metric_summ.output.scale,
        pair=rules.metric_summ.output.pair,
        stats=rules.topo_aggr.output,
    output: 'plt/eval/fig.pdf'
    shell:
        """
        python workflow/scripts/plt/eval/ranking_figure.py -i {input.metrics} -s {input.scale} -p {input.pair} -t {input.stats} -o {output}
        """
