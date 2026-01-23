localrules: annotate_eye

rule download_eye:
    threads: 4
    singularity: 'workflow/envs/gretabench.sif'
    input:
        img='workflow/envs/gretabench.sif',
        gid=rules.gen_gid_ensmbl.output.hg38,
    output:
        annot=temp(local('dts/hg38/eye/annot.csv')),
        frag=temp(local(expand('dts/hg38/eye/{sample}.frags.tsv.gz', sample=config['dts']['eye']['samples']))),
        rna=temp(local('dts/hg38/eye/rna.h5ad')),
    params:
        url_tar=config['dts']['eye']['url']['tar'],
        url_ann=config['dts']['eye']['url']['ann'],
        samples=config['dts']['eye']['samples'],
    shell:
        """
        path_eye=$(dirname {output.annot})
        path_tmp=$path_eye/tmp.tar
        wget --no-verbose "{params.url_ann}" -O {output.annot}
        wget --no-verbose "{params.url_tar}" -O $path_tmp
        tar -xvf $path_tmp -C $path_eye
        rm $path_eye/*.hic
        rm $path_eye/*.bedpe.gz
        python workflow/scripts/dts/eye/split.py \
        -a {output.annot} \
        -s {params.samples}
        samples=({params.samples})
        for sample in $path_eye/*.frags.tsv; do
            base_name=$(basename "$sample" .frags.tsv)
            # Check if sample is in configured samples list
            keep=false
            for s in "${{samples[@]}}"; do
                if [[ "$base_name" == "$s" ]]; then
                    keep=true
                    break
                fi
            done
            if [[ "$keep" == true ]]; then
                echo $sample 1>&2
                gzip $sample
            else
                rm $sample
            fi
        done
        python workflow/scripts/dts/eye/prcrna.py \
        {input.gid} \
        {output.annot} \
        {output.rna}
        rm $path_eye/GSM5866081_*
        rm $path_eye/GSM5866073_*
        rm $path_tmp
        """


rule callpeaks_eye:
    threads: 12
    singularity: 'workflow/envs/gretabench.sif'
    input:
        frags=rules.download_eye.output.frag,
        annot=rules.download_eye.output.annot,
    output: peaks=temp(local('dts/hg38/eye/peaks.h5ad'))
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


rule annotate_eye:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        rna=rules.download_eye.output.rna,
        peaks=rules.callpeaks_eye.output.peaks,
        annot=rules.download_eye.output.annot,
    output: out='dts/hg38/eye/annotated.h5mu'
    shell:
        """
        python workflow/scripts/dts/eye/eye.py \
        -a {input.rna} \
        -b {input.peaks} \
        -c {input.annot} \
        -d {output.out}
        """
