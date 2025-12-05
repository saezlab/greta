localrules: fig_dbs, fig_dbs_qc


rule fig_dbs:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        sts='anl/dbs/stats.csv',
        ovc='anl/dbs/ocoef.csv',
    output: 'plt/dbs/fig.pdf'
    shell:
        """
        python workflow/scripts/plt/dbs/stats.py {input.sts} {input.ovc} {output}
        """

rule fig_dbs_qc:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        knocktf=rules.prt_knocktf.output.dir,
        hpa=rules.tfm_hpa.output,
    output: 'plt/dbs/qc.pdf'
    shell:
        """
        python workflow/scripts/plt/dbs/qc.py {input.knocktf} {input.hpa} {output}
        """
