# heartatlas
rule download_fragments:
    output:
        tar=temp(local('datasets/heartatlas/fragments.tar')),
        frag=local(expand('datasets/heartatlas/{batchID}_atac_fragments.tsv.gz', batchID=config['datasets']['heartatlas']['batchIDs']))
    params:
        tar=config['datasets']['heartatlas']['url']['tar']
    shell:
        """
        data_path=$(dirname "{output.tar}")
        echo "Data path: $data_path"

        echo "Downloading tar file from {params.tar} to {output.tar}"
        wget '{params.tar}' -O '{output.tar}'

        echo "Extracting tar file to $data_path"
        tar -xvf '{output.tar}' -C "$data_path"

        echo "Removing .tbi files"
        rm "$data_path"/*.tbi
        """

rule download_anndata:
    output:
        anndata=local('datasets/heartatlas/multiome_raw.h5ad')
    params:
        anndata=config['datasets']['heartatlas']['url']['anndata']
    shell:
        """
        echo "Downloading anndata file from {params.anndata} to {output.anndata}"
        wget '{params.anndata}' -O '{output.anndata}'
        """

rule prcannot_heartatlas:
    input:
        h5ad='datasets/heartatlas/multiome_raw.h5ad'
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        annot='datasets/heartatlas/annot.csv'
    shell:
        """
        python workflow/scripts/datasets/heartatlas/heart_annot.py \
        -i {input.h5ad} \
        -a {output.annot}
        """

rule callpeaks_heartatlas:
    threads: 32
    input:
        frags=expand('datasets/heartatlas/{batchID}_atac_fragments.tsv.gz', batchID=config['datasets']['heartatlas']['batchIDs']),
        annot='datasets/heartatlas/annot.csv',
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        tmp=temp(directory(local('datasets/heartatlas/tmp_peaks'))),
        peaks=local('datasets/heartatlas/peaks.h5ad')
    resources:
        mem_mb=64000,
        runtime=2160,
    shell:
        """
        python workflow/scripts/datasets/callpeaks.py \
        -f {input.frags} \
        -a {input.annot} \
        -t {output.tmp} \
        -o {output.peaks}
        """

rule annotate_heartatlas:
    input:
        path_h5ad='datasets/heartatlas/multiome_raw.h5ad',
        path_peaks='datasets/heartatlas/peaks.h5ad',
        path_annot='datasets/heartatlas/annot.csv',
        path_geneids='gdata/geneids/',
    output:
        'datasets/heartatlas/annotated.h5mu'
    params:
        organism=config['datasets']['reprofibro']['organism']
    shell:
        """
        python workflow/scripts/datasets/reprofibro/reprofibro.py \
        -a {input.path_h5ad} \
        -b {input.path_annot} \
        -c {input.path_geneids} \
        -d {params.organism} \
        -e {input.path_peaks} \
        -f {output}
        """
