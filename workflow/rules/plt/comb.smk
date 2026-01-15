#localrules: fig_comb


rule fig_comb:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        mdta='dts/hg38/pitupair/cases/all/mdata.h5mu',
        qc='anl/dts/hg38.pitupair.all.qc.csv',
        nc='anl/dts/hg38.pitupair.all.nc.csv',
        sims='anl/topo/hg38.pitupair.all.sims_mult.csv',
        stat='anl/topo/hg38.pitupair.all.stats_mult.csv',
        fvsd='anl/topo/hg38.pitupair.all.fvsd.csv',
        stab='anl/stab/hg38.pitupair.all.ovsd.csv',
    output: 'plt/comb/fig.pdf'
    shell:
        """
        python workflow/scripts/plt/comb/sims.py {input.mdta} \
        {input.nc} {input.qc} {input.sims} {input.stat} {input.fvsd} {input.stab} {output}
        """
