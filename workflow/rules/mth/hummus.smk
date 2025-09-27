rule mdl_o_hummus:
    threads: 1
    singularity: 'workflow/envs/hummus.sif'
    input:
        mdata=rules.extract_case.output.mdata,
       # gen=rules.gen_genome_hummus.output
    output:
        rna_layer = 'dts/{dat}/cases/{case}/runs/tmp_o_hummus/rna_layer.tsv',
        atac_layer = 'dts/{dat}/cases/{case}/runs/tmp_o_hummus/atac_layer.tsv',
        out='dts/{dat}/cases/{case}/runs/o_hummus.o_hummus.o_hummus.o_hummus.mdl.csv'
    params:
        organism=lambda w: config['dts'][w.dat]['organism'],
        num_workers=24
    resources:
#        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'] * 2,
    shell:
        """
        # Extract annot
        path_tmp=$(dirname {output.out})/tmp_o_hummus
        mkdir -p $path_tmp
        # Create omics layers
        python workflow/scripts/mth/hummus/create_layers.py \
        -f {input.mdata} \
        --rna_layer {output.rna_layer} \
        --atac_layer {output.atac_layer} \
        --organism {params.organism} \
        -d $path_tmp \
        -c {params.num_workers}
        Rscript workflow/scripts/mth/hummus/src.R \
        {input.mdata} \
        {output.rna_layer} \
        {output.atac_layer} \
        $path_tmp \
        1 \
        {output.out}
        """
