rule mdl_o_scmtni:
    threads: 32
    singularity: 'workflow/envs/scmtni.sif'
    input:
        mdata=rules.extract_case.output.mdata,
        chainfiles=rules.download_liftover_chains.output.hg38ToHg19,
        motifs=rules.gen_motif_scmtni.output.motifs_dir,
        promoters=rules.gen_motif_scmtni.output.promoters_dir,
    output:
        out='dts/{dat}/cases/{case}/runs/o_scmtni.o_scmtni.o_scmtni.o_scmtni.mdl.csv'
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'] * 2,
    shell:
        """
        python workflow/scripts/mth/scmtni/src.py \
        -a {input.mdata} \
        -b {output.out} \
        -c {input.chainfiles} \
        -d {input.motifs} \
        -e {input.promoters}
        """
