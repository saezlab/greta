localrules: prcannot_heart


rule download_fragments_heart:
    threads: 7
    singularity: 'workflow/envs/figr.sif'
    input: 'workflow/envs/figr.sif'
    output:
        tar=temp(local('dts/hg38/heart/fragments.tar')),
        frag=expand('dts/hg38/heart/{sample}.frags.tsv.gz', sample=config['dts']['heart']['samples']),
        tbis=expand('dts/hg38/heart/{sample}.frags.tsv.gz.tbi', sample=config['dts']['heart']['samples'])
    params:
        tar=config['dts']['heart']['url']['tar']
    shell:
        """
        data_path=$(dirname "{output.tar}")
        wget --no-verbose '{params.tar}' -O '{output.tar}'
        tar -xvf '{output.tar}' -C "$data_path"
        rm "$data_path"/*.tbi
        rm {output.tar}
        for file in $data_path/*_atac_fragments.tsv.gz; do
            base_name=$(basename "$file" _atac_fragments.tsv.gz);
            new_file="${{base_name#*_}}.frags.tsv.gz";
            mv $file $data_path/$new_file
        done && \
        ls $data_path/*.frags.tsv.gz | xargs -n 1 -P {threads} bash workflow/scripts/dts/format_frags.sh
        """


rule download_anndata_heart:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: 'workflow/envs/gretabench.sif'
    output:
        adata=temp(local('dts/hg38/heart/multiome_raw.h5ad')),
        annot=temp(local('dts/hg38/heart/atac.h5ad'))
    params:
        adata=config['dts']['heart']['url']['anndata'],
        annot=config['dts']['heart']['url']['annot']
    shell:
        """
        wget --no-verbose '{params.adata}' -O '{output.adata}'
        wget --no-verbose '{params.annot}' -O '{output.annot}'
        """


rule prcannot_heart:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: rules.download_anndata_heart.output.annot,
    output: temp(local('dts/hg38/heart/annot.csv'))
    shell:
        """
        python workflow/scripts/dts/heart/heart_annot.py \
        -i {input} \
        -o {output}
        """


rule callpeaks_heart:
    threads: 8
    singularity: 'workflow/envs/gretabench.sif'
    input:
        frags=rules.download_fragments_heart.output.frag,
        annot=rules.prcannot_heart.output,
    output: peaks=temp(local('dts/hg38/heart/peaks.h5ad'))
    resources:
        mem_mb=128000,
        runtime=2160,
    shell:
        """
        python workflow/scripts/dts/callpeaks.py \
        -f {input.frags} \
        -a {input.annot} \
        -t $TMPDIR \
        -n {threads} \
        -o {output.peaks}
        """


rule annotate_heart:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        path_h5ad=rules.download_anndata_heart.output.adata,
        path_peaks=rules.callpeaks_heart.output.peaks,
        path_annot=rules.prcannot_heart.output,
        gid=rules.gen_gid_ensmbl.output
    output: out='dts/hg38/heart/annotated.h5mu'
    shell:
        """
        python workflow/scripts/dts/heart/heart.py \
        -a {input.path_h5ad} \
        -b {input.path_peaks} \
        -c {input.path_annot} \
        -e {input.gid} \
        -f {output}
        """
