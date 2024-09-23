rule download_tfsscenic:
    output:
        'gdata/tfs/scenic.txt'
    params:
        url='https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt'
    resources:
        mem_mb=2000
    shell:
        "wget '{params.url}' -O {output}"


rule grn_scenic:
    threads: 32
    input:
        data='datasets/{dataset}/cases/{case}/mdata.h5mu',
        tf='gdata/tfs/scenic.txt',
    	proms='/mnt/sds-hd/sd22b002/projects/GRETA/greta_resources/database/hg38/cre/promoters/promoters.bed'
    singularity:
    	'workflow/envs/scenicplus.sif'
    output:
        adj=temp(local('datasets/{dataset}/cases/{case}/runs/adj_tmp.tsv')),
    	t=temp(local('datasets/{dataset}/cases/{case}/runs/scenic_tmp.loom')),
    	grn='datasets/{dataset}/cases/{case}/runs/scenic.scenic.scenic.scenic.grn.csv'
    resources:
        mem_mb=64000
    shell:
        """
        # Step 1: Create Loom file
        python workflow/scripts/methods/scenic/loom.py \
        -i {input.data} \
        -o {output.t}

        # Step 2: Run pyscenic GRN
        pyscenic grn {output.t} {input.tf} -o {output.adj} --num_workers 30

    	# Step 3: Process GRN
    	python workflow/scripts/methods/scenic/process_grn.py \
    	-i {input.data} \
    	-o {output.grn} \
    	-p {input.proms} \
    	-g {output.adj} \
        """

