rule download_brain:
    threads: 34
    singularity: 'workflow/envs/figr.sif'
    output:
        tar=temp(local('dts/brain/GSE193688.tar')),
        annot=temp(local('dts/brain/raw_annot.csv')),
        frags=expand('dts/brain/{sample}.frags.tsv.gz', sample=config['dts']['brain']['samples']),
        tbis=expand('dts/brain/{sample}.frags.tsv.gz.tbi', sample=config['dts']['brain']['samples']),
        gex=temp(local(expand('dts/brain/{sample}_filtered_feature_bc_matrix.h5', sample=config['dts']['brain']['samples']))),
    params:
        full_dataset=config['dts']['brain']['url']['full_dataset'],
        annot=config['dts']['brain']['url']['annot'],
    shell:
        """
        data_path=$(dirname {output.tar})
        echo "Downloading tar file"
        wget --no-verbose '{params.full_dataset}' -O '{output.tar}'
        wget --no-verbose '{params.annot}' -O '{output.annot}'
        tar -xvf '{output.tar}' -C $data_path
        for file in $data_path/*_atac_fragments.tsv.gz; do
            base_name=$(basename "$file" _atac_fragments.tsv.gz);
            new_file="${{base_name#*_}}.frags.tsv.gz";
            mv $file $data_path/$new_file
        done && \
        ls $data_path/*.frags.tsv.gz | xargs -n 1 -P {threads} bash workflow/scripts/dts/format_frags.sh
        rm $data_path/*peaks.bed.gz
        (cd $data_path && for x in GSM*; do    mv $x `echo $x | cut -c 12-`; done)
        """


rule prc_annot:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: rules.download_brain.output.annot,
    output: annot=temp(local('dts/brain/annot.csv')),
    params: samples=config['dts']['brain']['samples'],
    shell:
        """
        python workflow/scripts/dts/brain/prc_annot.py \
        -a {input} \
        -b {params.samples} \
        -c {output.annot}
        """


rule callpeaks_brain:
    threads: 32
    singularity: 'workflow/envs/gretabench.sif'
    input:
        frags=rules.download_brain.output.frags,
        annot=rules.prc_annot.output.annot,
    output: peaks=temp(local('dts/brain/peaks.h5ad'))
    resources:
        mem_mb=110000,
        runtime=2160,
    shell:
        """
        python workflow/scripts/dts/callpeaks.py \
        -f {input.frags} \
        -a {input.annot} \
        -t '/tmp/brain/' \
        -n {threads} \
        -o {output.peaks}
        """


rule annotate_brain:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        path_gex=rules.download_brain.output.gex,
        path_peaks=rules.callpeaks_brain.output.peaks,
        path_annot=rules.prc_annot.output.annot,
        gid=rules.gen_gid_ensmbl.output,
    output: out='dts/brain/annotated.h5mu'
    shell:
        """
        python workflow/scripts/dts/brain/brain.py \
        -a {input.path_gex} \
        -b {input.path_peaks} \
        -c {input.path_annot} \
        -d {input.gid} \
        -f {output.out}
        """
