localrules: plt_fake


rule plt_fake:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        knn='anl/pair/pitu.all.fake_knn.csv',
        ctp='anl/pair/pitu.all.fake_prp.csv',
        cor='anl/pair/pitu.all.fake_cor.csv',
    output: 'plt/fig3/fake.pdf'
    shell:
        """
        python workflow/scripts/plt/fig3/fake.py {input.knn} {input.ctp} {input.cor} {output}
        """
