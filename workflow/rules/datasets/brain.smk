rule download_brain:
    threads: 34
    singularity: 'workflow/envs/figr.sif'
    output:
        tar=temp(local('datasets/brain/GSE193688.tar')),
        annot=temp(local('datasets/brain/raw_annot.csv')),
        frags=expand('datasets/brain/{sample}.frags.tsv.gz', sample=config['datasets']['brain']['samples']),
        tbis=expand('datasets/brain/{sample}.frags.tsv.gz.tbi', sample=config['datasets']['brain']['samples']),
        gex=temp(local(expand('datasets/brain/{sample}_filtered_feature_bc_matrix.h5', sample=config['datasets']['brain']['samples']))),
    params:
        full_dataset=config['datasets']['brain']['url']['full_dataset'],
        annot=config['datasets']['brain']['url']['annot'],
    shell:
        """
        data_path=$(dirname {output.tar})
        echo "Downloading tar file"
        wget --no-verbose '{params.full_dataset}' -O '{output.tar}'
        wget --no-verbose '{params.annot}' -O '{output.annot}'
        tar -xvf '{input.tar}' -C $data_path
        for file in $data_path/*_atac_fragments.tsv.gz; do
            base_name=$(basename "$file" _atac_fragments.tsv.gz);
            new_file="${{base_name#*_}}.frags.tsv.gz";
            mv $file $data_path/$new_file
        done && \
        ls $data_path/*.frags.tsv.gz | xargs -n 1 -P {threads} bash workflow/scripts/datasets/format_frags.sh
        rm $data_path/*peaks.bed.gz
        (cd $data_path && for x in GSM*; do    mv $x `echo $x | cut -c 12-`; done)
        """


rule prc_annot:
    threads: 1
    input:
        raw_annot=rules.download_brain.output.annot,
    output:
        annot=temp(local('datasets/brain/annot.csv')),
    singularity:
        'workflow/envs/gretabench.sif'
    params:
        samples=config['datasets']['brain']['samples'],
    shell:
        """
        python workflow/scripts/datasets/brain/prc_annot.py \
        -a {input.raw_annot} \
        -b {params.samples} \
        -c {output.annot}
        """


rule callpeaks_brain:
    threads: 32
    input:
        frags=rules.download_brain.output.frags,
        annot=rules.prc_annot.output.annot,
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
        -n {threads} \
        -o {output.peaks}
        """


rule annotate_brain:
    threads: 1
    input:
        path_gex=rules.download_brain.output.gex,
        path_peaks=rules.callpeaks_brain.output.peaks,
        path_annot=rules.prc_annot.output.annot,
        path_geneids=rules.download_geneids.output.dr,
    output:
        out='datasets/brain/annotated.h5mu'
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
        -f {output.out}
        """
