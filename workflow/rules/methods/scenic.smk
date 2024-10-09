localrules: download_tfsscenic, download_rankings


rule download_tfsscenic:
    output:
        out='gdata/tfs/scenic.txt'
    params:
        url='https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt'
    resources:
        mem_mb=2000
    shell:
        "wget '{params.url}' -O {output.out}"


rule download_rankings:
    output:
        small="aertslab/cistarget/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
    	big="aertslab/cistarget/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
    params:
        url_small="https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
    	url_big="https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
    shell:
        """
    	wget {params.url_big} -O {output.big}
    	wget {params.url_small} -O {output.small}
    	"""


rule mdl_scenic:
    threads: 16
    input:
        mdata=rules.extract_case.output.mdata,
        tf=rules.download_tfsscenic.output.out,
        proms='/mnt/sds-hd/sd22b002/projects/GRETA/greta_resources/database/hg38/cre/promoters/promoters.bed',
        ranking_small=rules.download_rankings.output.small,
        ranking_big=rules.download_rankings.output.big,
        motifs=rules.download_motifs.output.h
    singularity:
        'workflow/envs/scenicplus.sif'
    output:
        adj=temp(local('datasets/{dataset}/cases/{case}/runs/adj_tmp.tsv')),
        t=temp(local('datasets/{dataset}/cases/{case}/runs/scenic_tmp.loom')),
    	reg=temp(local('datasets/{dataset}/cases/{case}/runs/scenic_reg.csv')),
        out='datasets/{dataset}/cases/{case}/runs/scenic.scenic.scenic.scenic.mdl.csv'
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    shell:
        """
        # Step 1: Create Loom file
        python workflow/scripts/methods/scenic/loom.py \
        -i {input.mdata} \
        -o {output.t}
        echo "Created loom"

        # Step 2: Run pyscenic GRN
        pyscenic grn {output.t} {input.tf} -o {output.adj} --num_workers 16
        echo "Generated adj"

    	# Step 3: Run CTX
    	pyscenic ctx {output.adj} \
    	{input.ranking_small} \
    	{input.ranking_big} \
    	--annotations_fname {input.motifs} \
    	--expression_mtx_fname {output.t} \
    	--output {output.reg} \
    	--mask_dropouts \
    	--num_workers 16
        echo "Filtered TFs by motifs"

        # Step 4: Process GRN
        python workflow/scripts/methods/scenic/process_grn.py \
        -o {output.out} \
        -p {input.proms} \
        -g {output.adj} \
    	-r {output.reg}
        echo "Done"
        """
