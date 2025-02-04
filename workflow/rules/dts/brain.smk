localrules: prc_annot


rule download_brain:
    threads: 8
    singularity: 'workflow/envs/figr.sif'
    output:
        annot=temp(local('dts/brain/raw_annot.csv')),
        frags=expand('dts/brain/{sample}.frags.tsv.gz', sample=config['dts']['brain']['samples']),
        tbis=expand('dts/brain/{sample}.frags.tsv.gz.tbi', sample=config['dts']['brain']['samples']),
        gex=temp(local(expand('dts/brain/{sample}_filtered_feature_bc_matrix.h5', sample=config['dts']['brain']['samples']))),
    params:
        dtsts=config['dts']['brain']['url']['full_dataset'],
        annot=config['dts']['brain']['url']['annot'],
        samples=config['dts']['brain']['samples'],
    resources:
        mem=4000,
    shell:
        """
        allowed_samples=( {params.samples} )
        data_path=$(dirname {output.annot})
        curl -s {params.dtsts} | \
        grep -oP '(?<=acc=)GSM[0-9]+|(?<=<td valign="top">)[A-Za-z0-9]+' | \
        grep -v '^Illumina' | \
        paste -d'\t' - - | \
        while read gsm sample; do
          if [[ " ${{allowed_samples[@]}} " =~ " $sample " ]]; then
            base_path="${{gsm:0:7}}nnn"
            # Download atac frags
            url_frag="https://ftp.ncbi.nlm.nih.gov/geo/samples/${{base_path}}/${{gsm}}/suppl/${{gsm}}%5F${{sample}}%5Fatac%5Ffragments%2Etsv%2Egz"
            echo "$url_frag";
            wget -q "$url_frag" -O "$data_path/${{sample}}.frags.tsv.gz"
            # Download rna counts
            url_rna="https://www.ncbi.nlm.nih.gov/geo/download/?acc=${{gsm}}&format=file&file=${{gsm}}%5F${{sample}}%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
            echo "$url_rna";
            wget -q "$url_rna" -O "$data_path/${{sample}}_filtered_feature_bc_matrix.h5"
            echo "Downloaded ${{sample}}"
          fi
        done && \
        ls $data_path/*.frags.tsv.gz | xargs -n 1 -P {threads} bash workflow/scripts/dts/format_frags.sh && \
        wget --no-verbose '{params.annot}' -O '{output.annot}'
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
    threads: 8
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
        -t $TMPDIR \
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
