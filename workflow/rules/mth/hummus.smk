rule mdl_o_hummus:
    threads: 16
    singularity: 'workflow/envs/hummus.sif'
    input:
        img='workflow/envs/hummus.sif',
        mdata=rules.extract_case.output.mdata,
    output:
        out='dts/{org}/{dat}/cases/{case}/runs/o_hummus.o_hummus.o_hummus.o_hummus.mdl.csv'
    params:
        organism=lambda w: config['dts'][w.dat]['organism'],
        ext=config['methods']['hummus']['ext'],
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'] * 2,
    shell:
        """
        # Extract annot
        path_tmp=$(dirname {output.out})/tmp_o_hummus
        mkdir -p $path_tmp
        # Create omics layers
        python workflow/scripts/mth/hummus/create_layers.py \
        -f {input.mdata} \
        --rna_layer $path_tmp/rna_layer.tsv \
        --atac_layer $path_tmp/atac_layer.tsv \
        --organism {params.organism} \
        -d $path_tmp \
        -w {params.ext} \
        -c {threads}
        Rscript workflow/scripts/mth/hummus/src.R \
        {input.mdata} \
        $path_tmp/rna_layer.tsv \
        $path_tmp/atac_layer.tsv \
        $path_tmp \
        1 \
        $path_tmp/intermediar_grn.csv \
        $path_tmp/binding_regions.csv

        python workflow/scripts/mth/hummus/add_enhancers.py \
        -b $path_tmp/binding_regions.csv \
        -t $path_tmp/multilayer/bipartite/atac_rna.tsv \
        -r $path_tmp/rna_layer.tsv \
        -a $path_tmp/atac_layer.tsv \
        -i $path_tmp/intermediar_grn.csv \
        -o {output.out} \
        -g $path_tmp/genome_annotations.csv
        rm -rf $path_tmp
        """
