#localrules: fig_comb


rule fig_comb:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        mdta='dts/pbmc10k/cases/all/mdata.h5mu',
        qc='anl/dts/pbmc10k.all.qc.csv',
        nc='anl/dts/pbmc10k.all.nc.csv',
        sims='anl/topo/pbmc10k.all.sims_mult.csv',
        stat='anl/topo/pbmc10k.all.stats_mult.csv',
        fvsd='anl/topo/pbmc10k.all.fvsd.csv',
        stab='anl/stab/pbmc10k.all.ovsd.csv',
    output: 'plt/comb/fig.pdf'
    shell:
        """
        python workflow/scripts/plt/comb/sims.py {input.mdta} \
        {input.nc} {input.qc} {input.sims} {input.stat} {input.fvsd} {input.stab} {output}
        """
