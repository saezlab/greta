localrules: plt_pair


rule plt_pair:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        pmd='dts/pitupair/cases/all/mdata.h5mu',
        nmd='dts/pitunpair/cases/all/mdata.h5mu',
        ral='anl/pair/pitu.all.real_corvals.csv',
        qc='anl/pair/pitu.all.qc.csv',
        nc='anl/pair/pitu.all.ncells.csv',
        oc='anl/pair/pitu.all.pvsn.csv',
    output: 'plt/fig2/pair.pdf'
    shell:
        """
        python workflow/scripts/plt/fig2/pair.py {input.pmd} {input.nmd} {input.ral} {input.qc} {input.nc} {input.oc} {output}
        """
