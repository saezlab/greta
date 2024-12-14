localrules: plt_figs


rule plt_figs:
    threads: 1
    input: ['plt/fig1/fig1.pdf', 'plt/fig2/fig2.pdf', 'plt/fig3/fig3.pdf', 'plt/fig4/fig4.pdf', 'plt/fig5/fig5.pdf', 'plt/fig6/fig6.pdf']
    output: 'plt/figs.txt'
    shell:
        """
        touch {output}
        echo 'Done'
        """