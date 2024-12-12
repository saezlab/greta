localrules: plt_fig6


rule plt_fig6:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: 'anl/metrics/summary/pbmc10k.all.csv'
    output: 'plt/fig6/fig6.pdf'
    shell:
        """
        python workflow/scripts/plt/fig6/eval.py {input} {output}
        """
