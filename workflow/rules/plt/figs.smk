localrules: plt_figs


rule plt_figs:
    threads: 1
    input: ['plt/stab/fig.pdf', 'plt/pair/fig.pdf', 'plt/comb/fig.pdf', 'plt/dbs/fig.pdf', 'plt/eval/fig.pdf']
    output: 'plt/figs.txt'
    shell:
        """
        touch {output}
        echo 'Done'
        """