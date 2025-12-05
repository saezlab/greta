localrules: prcannot_lung, annotate_lung


rule download_lung:
    threads: 4
    singularity: 'workflow/envs/figr.sif'
    input: 'workflow/envs/figr.sif'
    output:
        raw_annot=temp(local('dts/hg38/lung/raw_annot.csv')),
        frag=temp(local(expand('dts/hg38/lung/{sample}.frags.tsv.gz', sample=config['dts']['lung']['samples']))),
        gex=temp(local(expand('dts/hg38/lung/{sample}_matrix.h5', sample=config['dts']['lung']['samples']))),
    params:
        url_geo=config['dts']['lung']['url']['geo'],
        url_annot=config['dts']['lung']['url']['annot'],
        gex_geoids=config['dts']['lung']['gex_geoids'],
        atac_geoids=config['dts']['lung']['atac_geoids'],
        samples=config['dts']['lung']['samples'],
    shell:
        """
        path_lung=$(dirname {output.raw_annot})
        wget --no-verbose "{params.url_annot}" -O {output.raw_annot}
        for GEOID_SAMPLE in {params.gex_geoids}; do
            SAMPLE=$(echo "${{GEOID_SAMPLE#*%5F}}")
            GEOID=$(echo "${{GEOID_SAMPLE%%'%5F'*}}")
            URL_RNA=$(echo "{params.url_geo}/?acc=${{GEOID}}&format=file&file=${{GEOID_SAMPLE}}%5Fraw%5Ffeature%5Fbc%5Fmatrix.h5")
            echo $GEOID $SAMPLE
            wget --no-verbose "$URL_RNA" -O $path_lung/${{SAMPLE}}_matrix.h5
        done
        for GEOID_SAMPLE in {params.atac_geoids}; do
            SAMPLE=$(echo "${{GEOID_SAMPLE#*%5F}}")
            GEOID=$(echo "${{GEOID_SAMPLE%%'%5F'*}}")
            URL_ATAC=$(echo "{params.url_geo}/?acc=${{GEOID}}&format=file&file=${{GEOID_SAMPLE}}%5Fatac%5Ffragments.tsv.gz")
            echo $GEOID $SAMPLE
            wget --no-verbose "$URL_ATAC" -O $path_lung/${{SAMPLE}}.frags.tsv.gz
        done
        ls $path_lung/*.frags.tsv.gz | xargs -n 1 -P {threads} bash workflow/scripts/dts/format_frags.sh
        """


rule prcannot_lung:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        raw_annot=rules.download_lung.output.raw_annot,
        gex=rules.download_lung.output.gex,
        gid=rules.gen_gid_ensmbl.output.hg38,
    output:
        annot=temp(local('dts/hg38/lung/annot.csv')),
        rna=temp(local('dts/hg38/lung/rna.h5ad'))
    params:
        samples=config['dts']['lung']['samples'],
    shell:
        """
        python workflow/scripts/dts/lung/annot.py \
        -a {input.raw_annot} \
        -b {input.gid} \
        -c {params.samples} \
        -d {output.rna} \
        -e {output.annot}
        """


rule callpeaks_lung:
    threads: 16
    singularity: 'workflow/envs/gretabench.sif'
    input:
        frags=rules.download_lung.output.frag,
        annot=rules.prcannot_lung.output.annot,
    output: peaks=temp(local('dts/hg38/lung/peaks.h5ad'))
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

rule annotate_lung:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        gex=rules.prcannot_lung.output.rna,
        peaks=rules.callpeaks_lung.output.peaks,
        annot=rules.prcannot_lung.output.annot,
    output: out='dts/hg38/lung/annotated.h5mu'
    shell:
        """
        python workflow/scripts/dts/lung/lung.py \
        -a {input.gex} \
        -b {input.peaks} \
        -c {input.annot} \
        -d {output.out}
        """
