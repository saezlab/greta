rule mdl_o_crema:
    threads: 1
    singularity: 'workflow/envs/crema.sif'
    input:
        mdata=rules.extract_case.output.mdata,
        gen=rules.gen_genome_crema.output,
        tfs=rules.gen_motif_crema.output,
    output:
        out='dts/{dat}/cases/{case}/runs/o_crema.o_crema.o_crema.o_crema.mdl.csv'
    params:
        ext=config['methods']['crema']['ext'] // 2,
        site_extension=config['methods']['crema']['site_extension'],
        thr_fdr=config['methods']['crema']['thr_fdr'],
    shell:
        """
        # Extract annot
        path_tmp=$(dirname {output.out})/tmp_o_crema
        mkdir -p $path_tmp
        python -c "import mudata as mu; mu.read('{input.mdata}').obs.to_csv('$path_tmp/ann.csv');"
        # Create insertion file
        path_frags=($(dirname $(dirname $(dirname {input.mdata})))/*.frags.tsv.gz)
        python workflow/scripts/mth/crema/get_inserts.py \
        -f "${{path_frags[@]}}" \
        -g {input.gen} \
        -d $path_tmp
        bgzip $path_tmp/inserts.tsv
        tabix -p bed $path_tmp/inserts.tsv.gz
        Rscript workflow/scripts/mth/crema/src.R \
        {input.mdata} \
        $path_tmp/inserts.tsv.gz \
        {input.tfs} \ 
        {params.ext} \
        {params.site_extension} \
        {params.thr_fdr} \
        {output.out}
        """
