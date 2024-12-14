localrules: plt_fig4


rule plt_fig4:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        mdta='dts/pbmc10k/cases/all/mdata.h5mu',
        qc='anl/dts/pbmc10k.all.qc.csv',
        nc='anl/dts/pbmc10k.all.nc.csv',
        sims='anl/topo/pbmc10k.all.sims_mult.csv',
        stat='anl/topo/pbmc10k.all.stats_mult.csv',
        stab='anl/stab/pitupair.all.ovsd.csv',
    output: 'plt/fig4/fig4.pdf'
    shell:
        """
        python workflow/scripts/plt/fig4/sims.py {input.mdta} \
        {input.nc} {input.qc} {input.sims} {input.stat} {input.stab} {output}
        """
