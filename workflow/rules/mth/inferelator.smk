rule hvg_inf:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: rules.extract_case.output.mdata
    output: 'dts/{org}/{dat}/cases/{case}/runs/inferelator_hvg.txt.gz'
    shell:
        """
        python workflow/scripts/mth/inferelator/hvg.py {input} {output}
        """

rule mdl_o_inferelator:
    threads: 16
    singularity: 'workflow/envs/inferelator.sif'
    input:
        img='workflow/envs/inferelator.sif',
        mdata=rules.extract_case.output.mdata,
        gid=rules.gen_gid_ensmbl.output,
        fa=rules.gen_genome_inferelator.output.fa,
        gtf=rules.gen_genome_inferelator.output.gtf,
        meme=rules.gen_motif_inferelator.output.meme,
        hvg=rules.hvg_inf.output,
    output:
        out='dts/{org}/{dat}/cases/{case}/runs/o_inferelator.o_inferelator.o_inferelator.o_inferelator.mdl.csv'
    params:
        ext=config['methods']['inferelator']['ext'] // 2,
    resources:
        mem_mb=restart_mem,
        runtime=1440,
    shell:
        """
        path_tmp=$(dirname {output.out})/tmp_o_inferelator
        mkdir -p $path_tmp
        export TMPDIR=$path_tmp
        python workflow/scripts/mth/inferelator/prior.py \
        {input.mdata} \
        {input.hvg} \
        {input.fa} \
        {input.gtf} \
        {input.meme} \
        {params.ext} \
        $path_tmp
        python workflow/scripts/mth/inferelator/run.py \
        $path_tmp \
        {threads} \
        {output.out}
        rm -rf $path_tmp
        """
