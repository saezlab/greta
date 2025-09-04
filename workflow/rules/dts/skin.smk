localrules: download_skin, annotate_skin


rule download_skin:
    threads: 2
    singularity: 'workflow/envs/gretabench.sif'
    input:
        img='workflow/envs/gretabench.sif',
        gid=rules.gen_gid_ensmbl.output,
    output:
        ann=temp(local('dts/skin/annot.csv')),
        frag=temp(local(expand('dts/skin/{sample}.frags.tsv.gz', sample=config['dts']['skin']['samples']))),
        rna=temp(local('dts/skin/rna.h5ad')),
    params:
        url_tar=config['dts']['skin']['url']['tar'],
        url_ann=config['dts']['skin']['url']['ann'],
        samples=config['dts']['skin']['samples'],
    shell:
        """
        path_skin=$(dirname {output.ann})
        path_tmp=$path_skin/tmp.tar
        wget --no-verbose "{params.url_ann}" -O {output.ann}
        wget --no-verbose "{params.url_tar}" -O $path_tmp
        tar -xvf $path_tmp -C $path_skin
        rm $path_skin/*.tsv.gz.tbi.gz
        rm $path_skin/*_per_barcode_metrics.csv.gz
        for f in $path_skin/GSM*_atac_fragments.tsv.gz; do
            base=$(basename "$f")
            sample=${{base#*_}}
            sample=${{sample%%_atac_fragments.tsv.gz}}
            mv "$f" "$path_skin/${{sample}}.frags.tsv.gz"
        done
        python workflow/scripts/dts/skin/prcrna.py \
        {output.ann} \
        {input.gid}
        rm $path_skin/*_filtered_feature_bc_matrix.h5
        ls $path_skin/*.frags.tsv.gz | xargs -n 1 -P {threads} bash workflow/scripts/dts/format_frags.sh
        rm $path_tmp
        """


rule callpeaks_skin:
    threads: 8
    singularity: 'workflow/envs/gretabench.sif'
    input:
        frags=rules.download_skin.output.frag,
        annot=rules.download_skin.output.ann,
    output: peaks=temp(local('dts/skin/peaks.h5ad'))
    resources:
        mem_mb=64000,
        runtime=360,
    shell:
        """
        python workflow/scripts/dts/callpeaks.py \
        -f {input.frags} \
        -a {input.annot} \
        -t $TMPDIR \
        -n {threads} \
        -o {output.peaks}
        """


rule annotate_skin:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        rna=rules.download_skin.output.rna,
        peaks=rules.callpeaks_skin.output.peaks,
        annot=rules.download_skin.output.ann,
    output: out='dts/skin/annotated.h5mu'
    shell:
        """
        python workflow/scripts/dts/skin/skin.py \
        -a {input.rna} \
        -b {input.peaks} \
        -c {input.annot} \
        -d {output.out}
        """
