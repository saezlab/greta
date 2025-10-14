rule mdl_o_grnboost:
    threads: 16
    singularity: 'workflow/envs/scenicplus.sif'
    input:
        img='workflow/envs/scenicplus.sif',
        mdata=rules.extract_case.output.mdata,
        tf=rules.gen_tfs_scenic.output,
        proms=rules.cre_promoters.output,
    output:
        out='dts/{dat}/cases/{case}/runs/o_grnboost.o_grnboost.o_grnboost.o_grnboost.mdl.csv'
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'] * 2,
    shell:
        """
        # Create Loom file
        path_tmp=$(dirname {output.out})
        path_tmp=$path_tmp/tmp_grnboost
        path_gex=$path_tmp/gex.loom
        path_adj=$path_tmp/adj.tsv
        rm -rf $path_tmp
        mkdir -p $path_tmp
        python workflow/scripts/mth/scenic/loom.py \
        -i {input.mdata} \
        -o $path_gex
        echo "Created loom"

        # Run grnboost2 GRN
        arboreto_with_multiprocessing.py $path_gex {input.tf} -o $path_adj --num_workers {threads} --seed 42
        echo "Generated adj"

        # Process GRN
        python workflow/scripts/mth/scenic/process_grn.py \
        -g $path_adj \
        -p {input.proms} \
        -o {output.out}

        rm -rf $path_tmp
        """

rule mdl_o_scenic:
    threads: 16
    singularity: 'workflow/envs/scenicplus.sif'
    input:
        img='workflow/envs/scenicplus.sif',
        mdata=rules.extract_case.output.mdata,
        tf=rules.gen_tfs_scenic.output,
        proms=rules.cre_promoters.output,
        ranking_small=rules.gen_motif_scenic_rnk.output.sml,
        ranking_big=rules.gen_motif_scenic_rnk.output.big,
        motifs=rules.gen_motif_scenic.output
    output:
        out='dts/{dat}/cases/{case}/runs/o_scenic.o_scenic.o_scenic.o_scenic.mdl.csv'
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'] * 2,
    shell:
        """
        path_tmp=$(dirname {output.out})
        path_tmp=$path_tmp/tmp_scenic
        path_gex=$path_tmp/gex.loom
        path_adj=$path_tmp/adj.tsv
        path_reg=$path_tmp/reg.csv
        rm -rf $path_tmp
        mkdir -p $path_tmp

        python workflow/scripts/mth/scenic/loom.py \
        -i {input.mdata} \
        -o $path_gex
        echo "Created loom"

        arboreto_with_multiprocessing.py $path_gex {input.tf} -o $path_adj --num_workers {threads} --seed 42
        echo "Generated adj"

        # Step 3: Run CTX
        pyscenic ctx $path_adj \
        {input.ranking_small} \
        {input.ranking_big} \
        --annotations_fname {input.motifs} \
        --expression_mtx_fname $path_gex \
        --output $path_reg \
        --mask_dropouts \
        --num_workers {threads}
        echo "Filtered TFs by motifs"

        # Step 4: Process GRN
        python workflow/scripts/mth/scenic/process_grn.py \
        -g $path_adj \
        -p {input.proms} \
        -r $path_reg \
        -o {output.out}
        echo "Done"

        rm -rf $path_tmp
        """
