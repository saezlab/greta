
rule download_TF_list_scenic:
        output:
                'gdata/tfs/allTFs.txt'
        params:
                url='https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt'
        resources:
                mem_mb=2000
        shell:
                "wget '{params.url}' -O {output}"


rule grn_scenic:
    input:
        data='datasets/{dataset}/cases/{case}/mdata.h5mu',
        tf='gdata/tfs/allTFs.txt',
	proms='data/promoters.bed'   # we need to change location
    singularity:
	'workflow/envs/scenicplus.sif'
    output:
        adj=temp(local('datasets/{dataset}/cases/{case}/runs/adj_tmp.tsv')),     #is this okei?
	t=temp(local('datasets/{dataset}/cases/{case}/runs/scenic_tmp.loom')),
	grn='datasets/{dataset}/cases/{case}/runs/scenic.scenic.scenic.scenic.grn.csv'
    shell:
        """
        # Step 1: Create Loom file
        python workflow/scripts/methods/scenic/loom.py \
        -i {input.data} \
        -o {output.t}

        # Step 2: Run pyscenic GRN
        pyscenic grn {output.t} {input.tf} -o {output.adj}

	# Step 3: Process GRN
	python workflow/scripts/methods/scenic/process_grn.py \
	-i {input.data} \
	-o {output.grn} \
	-p {input.proms} \
	-g {output.adj} \
        """

