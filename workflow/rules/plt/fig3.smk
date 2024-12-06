localrules: plt_fig3


rule plt_fig3:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        knn='anl/pair/pitu.all.fake_knn.csv',
        ctp='anl/pair/pitu.all.fake_prp.csv',
        cor='anl/pair/pitu.all.fake_cor.csv',
        ocf='anl/pair/pitu.all.pvsf.csv',
    output: 'plt/fig3/fig3.pdf'
    shell:
        """
        python workflow/scripts/plt/fig3/fake.py {input.knn} {input.ctp} {input.cor} {input.ocf} {output}
        """
