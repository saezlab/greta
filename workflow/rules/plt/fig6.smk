localrules: plt_fig6


rule plt_fig6:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        smr='anl/metrics/summary/pbmc10k.all.csv',
        dct='anl/stab/dictys/pbmc10k.scores.csv',
    output: 'plt/fig6/fig6.pdf'
    shell:
        """
        python workflow/scripts/plt/fig6/eval.py {input} {output}
        """
