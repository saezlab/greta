rule mdl_o_scmtni:
    threads: 64
    singularity: 'workflow/envs/scmtni.sif'
    input:
        mdata=rules.extract_case.output.mdata,
        chainfiles=rules.download_liftover_chains.output.hg38ToHg19,
        motifs=rules.gen_motif_scmtni.output.motifs_dir,
        promoters=rules.gen_motif_scmtni.output.promoters_dir,
    output:
        out='dts/{dat}/cases/{case}/runs/o_scmtni.o_scmtni.o_scmtni.o_scmtni.mdl.csv'
    params:
        ext=config['methods']['scmtni']['ext'],
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'] * 2,
    shell:
        """
        path_out=$(dirname {output.out})
        path_out=$path_out/tmp_o_scmtni
        mkdir -p $path_out
        python workflow/scripts/mth/scmtni/src.py \
        -a {input.mdata} \
        -b $path_out \
        -c {output.out} \
        -d {input.chainfiles} \
        -e {input.motifs} \
        -f {input.promoters} \
        -g {params.ext} \
        -i {threads}
        """
