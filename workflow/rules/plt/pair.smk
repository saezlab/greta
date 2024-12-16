localrules: plt_npair, plt_fake, fig_pair


rule plt_npair:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        pmd='dts/pitupair/cases/all/mdata.h5mu',
        nmd='dts/pitunpair/cases/all/mdata.h5mu',
        ral='anl/pair/pitu.all.real_corvals.csv',
        qc='anl/pair/pitu.all.qc.csv',
        nc='anl/pair/pitu.all.ncells.csv',
        oc='anl/pair/pitu.all.pvsn.csv',
    output: 'plt/pair/npair.pdf'
    shell:
        """
        python workflow/scripts/plt/pair/pair.py {input.pmd} {input.nmd} {input.ral} {input.qc} {input.nc} {input.oc} {output}
        """


rule plt_fake:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        knn='anl/pair/pitu.all.fake_knn.csv',
        ctp='anl/pair/pitu.all.fake_prp.csv',
        cor='anl/pair/pitu.all.fake_cor.csv',
        ocf='anl/pair/pitu.all.pvsf.csv',
    output: 'plt/pair/fake.pdf'
    shell:
        """
        python workflow/scripts/plt/pair/fake.py {input.knn} {input.ctp} {input.cor} {input.ocf} {output}
        """


rule fig_pair:
    threads: 1
    input: ['plt/pair/npair.pdf', 'plt/pair/fake.pdf']
    output: 'plt/pair/fig.pdf'
    shell:
        """
        gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile={output} {input}
        """
