rule download_fragments_heart:
    threads: 1
    output:
        tar=temp(local('datasets/heartatlas/fragments.tar')),
        frag=expand('datasets/heartatlas/{sample}.frags.tsv.gz', sample=config['datasets']['heartatlas']['samples'])
    params:
        tar=config['datasets']['heartatlas']['url']['tar']
    shell:
        """
        data_path=$(dirname "{output.tar}")
        wget --no-verbose '{params.tar}' -O '{output.tar}'
        tar -xvf '{output.tar}' -C "$data_path"
        rm "$data_path"/*.tbi
        for file in $data_path/*_atac_fragments.tsv.gz; do
            new_file=$(echo "$file" | sed 's/_atac_fragments.tsv.gz/.frags.tsv.gz/')
            mv "$file" "$new_file"
            bash workflow/scripts/datasets/format_frags.sh $new_file
        done
        """


rule download_anndata_heart:
    threads: 1
    output:
        adata=temp(local('datasets/heartatlas/multiome_raw.h5ad')),
        annot=temp(local('datasets/heartatlas/atac.h5ad'))
    params:
        adata=config['datasets']['heartatlas']['url']['anndata'],
        annot=config['datasets']['heartatlas']['url']['annotation']
    shell:
        """
        wget --no-verbose '{params.adata}' -O '{output.adata}'
        wget --no-verbose '{params.annot}' -O '{output.annot}'
        """


rule prcannot_heartatlas:
    threads: 1
    singularity:
        'workflow/envs/gretabench.sif'
    input:
        h5ad=rules.download_anndata_heart.output.adata,
        atac=rules.download_anndata_heart.output.annot,
    output:
        annot=temp(local('datasets/heartatlas/annot.csv'))
    shell:
        """
        python workflow/scripts/datasets/heartatlas/heart_annot.py \
        -i {input.h5ad} \
        -a {input.atac} \
        -o {output.annot}
        """


rule callpeaks_heartatlas:
    threads: 32
    resources:
        mem_mb=512000,
        runtime=2160,
    singularity:
        'workflow/envs/gretabench.sif'
    input:
        frags=rules.download_fragments_heart.output.frag,
        annot=rules.prcannot_heartatlas.output.annot,
    output:
        tmp=temp(directory(local('datasets/heartatlas/tmp_peaks'))),
        peaks=temp(local('datasets/heartatlas/peaks.h5ad'))
    shell:
        """
        python workflow/scripts/datasets/callpeaks.py \
        -f {input.frags} \
        -a {input.annot} \
        -t {output.tmp} \
        -n {threads} \
        -o {output.peaks}
        """


rule annotate_heartatlas:
    threads: 1
    singularity:
        'workflow/envs/gretabench.sif'
    input:
        path_h5ad=rules.download_anndata_heart.output.adata,
        path_peaks=rules.callpeaks_heartatlas.output.peaks,
        path_annot=rules.prcannot_heartatlas.output.annot,
        path_geneids=rules.download_geneids.output.dr,
    output:
        out='datasets/heartatlas/annotated.h5mu'
    params:
        organism=config['datasets']['heartatlas']['organism']
    shell:
        """
        python workflow/scripts/datasets/heartatlas/heartatlas.py \
        -a {input.path_h5ad} \
        -b {input.path_annot} \
        -c {input.path_geneids} \
        -d {params.organism} \
        -e {input.path_peaks} \
        -f {output}
        """
