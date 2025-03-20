rule mdl_scenic:
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
        adj=temp(local('dts/{dat}/cases/{case}/runs/adj_tmp.tsv')),
        t=temp(local('dts/{dat}/cases/{case}/runs/scenic_tmp.loom')),
        reg=temp(local('dts/{dat}/cases/{case}/runs/scenic_reg.csv')),
        out='dts/{dat}/cases/{case}/runs/scenic.scenic.scenic.scenic.mdl.csv'
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'] * 2,
    shell:
        """
        # Step 1: Create Loom file
        python workflow/scripts/mth/scenic/loom.py \
        -i {input.mdata} \
        -o {output.t}
        echo "Created loom"

        # Step 2: Run pyscenic GRN
        arboreto_with_multiprocessing.py {output.t} {input.tf} -o {output.adj} --num_workers {threads} --seed 42
        echo "Generated adj"

        # Step 3: Run CTX
        pyscenic ctx {output.adj} \
        {input.ranking_small} \
        {input.ranking_big} \
        --annotations_fname {input.motifs} \
        --expression_mtx_fname {output.t} \
        --output {output.reg} \
        --mask_dropouts \
        --num_workers {threads}
        echo "Filtered TFs by motifs"

        # Step 4: Process GRN
        python workflow/scripts/mth/scenic/process_grn.py \
        -o {output.out} \
        -p {input.proms} \
        -g {output.adj} \
        -r {output.reg}
        echo "Done"
        """
