localrules: dts_qcstats


rule dts_qcstats:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: rules.extract_case.output.mdata
    output:
        qc='anl/dts/{dat}.{case}.qc.csv',
        nc='anl/dts/{dat}.{case}.nc.csv',
    shell:
        """
        python workflow/scripts/anl/dts/qcstats.py \
        {input} {output.qc} {output.nc}
        """
