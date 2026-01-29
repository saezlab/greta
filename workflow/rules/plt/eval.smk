localrules: fig_eval


rule fig_eval:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: rules.metric_summ.output.metrics
    output: 'plt/eval/fig.pdf'
    shell:
        """
        python workflow/scripts/plt/eval/ranking_figure.py -i {input} -o {output}
        """
