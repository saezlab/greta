localrules: fig_eval


rule fig_eval:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        smr='anl/metrics/summary/pbmc10k.all.csv',
        dct='anl/stab/unsmthds/pbmc10k.scores.csv',
    output: 'plt/eval/fig.pdf'
    shell:
        """
        python workflow/scripts/plt/eval/eval.py {input} {output}
        """
