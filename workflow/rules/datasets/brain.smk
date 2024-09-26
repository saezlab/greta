
SAMPLES = config['datasets']['brain']['samples']

rule download_brain:
    output:
        tar=local('datasets/brain/GSE193688.tar'), 
        annot=temp(local('datasets/brain/raw_annot.csv'))

    params:
        full_dataset=config['datasets']['brain']['url']['full_dataset'],
        annot=config['datasets']['brain']['url']['annot']


    shell:
        """
        data_path=$(dirname {output.tar})
 
        echo "Downloading tar file"
        wget '{params.full_dataset}' -O '{output.tar}'
        wget '{params.annot}' -O '{output.annot}'
    
        """

rule extract_files:
    input:
        tar=local('datasets/brain/GSE193688.tar')
    output:
        frag=temp(local(expand('datasets/brain/{sample}.frags.tsv.gz', sample=SAMPLES))),
        multi=temp(local(expand('datasets/brain/{sample}_filtered_feature_bc_matrix.h5', sample=SAMPLES)))
    shell:
        """
        data_path=$(dirname {input.tar})

        echo "Extracting files and removing archive"
        tar -xvf '{input.tar}' -C $data_path
        rm '{input.tar}'

        echo "Removing peak files"
        rm $data_path/*peaks.bed.gz

        echo "Rename files"
        (cd $data_path && for x in GSM*; do    mv $x `echo $x | cut -c 12-`; done)

        """


rule prc_annot:
    input:
        raw_annot=local('datasets/brain/raw_annot.csv'),
        multi=local(expand('datasets/brain/{sample}_filtered_feature_bc_matrix.h5', sample=SAMPLES))
    output:
        annot=local('datasets/brain/annot.csv')

    params:
        sample_dir=local('datasets/brain/')

    singularity:
        'workflow/envs/gretabench.sif'
        
    shell:
        """
        python workflow/scripts/datasets/brain/prc_annot.py \
        -a {input.raw_annot} \
        -b {output.annot} \
        -c {params.sample_dir} 
        """


rule callpeaks_brain:
    threads: 32
    input:
        frags=expand('datasets/brain/{sample}.frags.tsv.gz', sample=SAMPLES),
        annot='datasets/brain/annot.csv'
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        tmp=temp(directory(local('datasets/brain/tmp_peaks'))),
        peaks=temp(local('datasets/brain/peaks.h5ad'))
    resources:
        mem_mb=110000,
        runtime=2160,
    shell:
        """
        python workflow/scripts/datasets/callpeaks.py \
        -f {input.frags} \
        -a {input.annot} \
        -t {output.tmp} \
        -o {output.peaks}
        """


rule annotate_brain:
    input:
        path_gex=expand('datasets/brain/{sample}_filtered_feature_bc_matrix.h5', sample=SAMPLES),
        path_peaks='datasets/brain/peaks.h5ad',
        path_annot='datasets/brain/annot.csv',
        path_geneids='gdata/geneids/'
    output:
        'datasets/brain/annotated.h5mu'
    singularity:
        'workflow/envs/gretabench.sif'
    params:
        organism=config['datasets']['brain']['organism']
    shell:
        """
        python workflow/scripts/datasets/brain/brain.py \
        -a {input.path_gex} \
        -b {input.path_peaks} \
        -c {input.path_annot} \
        -d {input.path_geneids} \
        -e {params.organism} \
        -f {output}
        """
