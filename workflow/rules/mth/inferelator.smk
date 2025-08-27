rule mdl_o_inferelator:
    threads: 16
    singularity: 'workflow/envs/inferelator.sif'
    input:
        mdata=rules.extract_case.output.mdata,
        gid=rules.gen_gid_ensmbl.output,
        fa=rules.gen_genome_inferelator.output.fa,
        gtf=rules.gen_genome_inferelator.output.gtf,
        meme=rules.gen_motif_inferelator.output.meme,
    output:
        tmp=temp(local(directory('dts/{dat}/cases/{case}/runs/tmp_o_inferelator/'))),
        out='dts/{dat}/cases/{case}/runs/o_inferelator.o_inferelator.o_inferelator.o_inferelator.mdl.csv'
    params:
        ext=config['methods']['inferelator']['ext'] // 2,
    shell:
        """
        export TMPDIR={output.tmp}
        python workflow/scripts/mth/inferelator/prior.py \
        {input.mdata} \
        {input.fa} \
        {input.gtf} \
        {input.meme} \
        {params.ext} \
        {output.tmp}

        python workflow/scripts/mth/inferelator/run.py \
        {output.tmp} \
        {threads} \
        {output.out}
        """
