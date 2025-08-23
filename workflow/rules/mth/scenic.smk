rule mdl_grnboost:
    threads: 16
    singularity: 'workflow/envs/scenicplus.sif'
    input:
        img='workflow/envs/scenicplus.sif',
        mdata=rules.extract_case.output.mdata,
        tf=rules.gen_tfs_scenic.output,
        proms=rules.cre_promoters.output,
    output:
        gex=temp(local('dts/{dat}/cases/{case}/runs/tmp_scenic/gex.loom')),
        adj=temp(local('dts/{dat}/cases/{case}/runs/tmp_scenic/adj.tsv')),
        out='dts/{dat}/cases/{case}/runs/grnboost.grnboost.grnboost.grnboost.mdl.csv'
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'] * 2,
    shell:
        """
        # Create Loom file
        python workflow/scripts/mth/scenic/loom.py \
        -i {input.mdata} \
        -o {output.gex}
        echo "Created loom"

        # Run grnboost2 GRN
        arboreto_with_multiprocessing.py {output.gex} {input.tf} -o {output.adj} --num_workers {threads} --seed 42
        echo "Generated adj"

        # Process GRN
        python workflow/scripts/mth/scenic/process_grn.py \
        -g {output.adj} \
        -p {input.proms} \
        -o {output.out}
        """

rule mdl_scenic:
    threads: 16
    singularity: 'workflow/envs/scenicplus.sif'
    input:
        img='workflow/envs/scenicplus.sif',
        gex=rules.mdl_grnboost.output.gex,
        adj=rules.mdl_grnboost.output.adj,
        proms=rules.cre_promoters.output,
        ranking_small=rules.gen_motif_scenic_rnk.output.sml,
        ranking_big=rules.gen_motif_scenic_rnk.output.big,
        motifs=rules.gen_motif_scenic.output
    output:
        reg=temp(local('dts/{dat}/cases/{case}/runs/tmp_scenic/reg.csv')),
        out='dts/{dat}/cases/{case}/runs/scenic.scenic.scenic.scenic.mdl.csv'
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'] * 2,
    shell:
        """
        # Step 3: Run CTX
        pyscenic ctx {input.adj} \
        {input.ranking_small} \
        {input.ranking_big} \
        --annotations_fname {input.motifs} \
        --expression_mtx_fname {input.gex} \
        --output {output.reg} \
        --mask_dropouts \
        --num_workers {threads}
        echo "Filtered TFs by motifs"

        # Step 4: Process GRN
        python workflow/scripts/mth/scenic/process_grn.py \
        -g {input.adj} \
        -p {input.proms} \
        -r {output.reg} \
        -o {output.out}
        echo "Done"
        """
