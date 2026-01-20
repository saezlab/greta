localrules: prcannot_breast, annotate_breast


rule download_fragments_breast:
    threads: 2
    singularity: 'workflow/envs/figr.sif'
    input: 'workflow/envs/figr.sif'
    output:
        frag=expand('dts/hg38/breast/{sample}.frags.tsv.gz', sample=config['dts']['breast']['samples'])
    params:
        tar=config['dts']['breast']['url']['tar'],
        samples=config['dts']['breast']['samples'],
    resources:
        mem_mb=8000,
    shell:
        """
        paths=({output.frag})
        path_breast=$(dirname "${{paths[0]}}")
        path_tmp=$path_breast/tmp.tsv.bgz
        wget --no-verbose '{params.tar}' -O $path_tmp
        # Process each pool in parallel (one core per pool)
        parallel -j{threads} workflow/scripts/dts/breast/write_sample.sh {{}} "$path_breast" "$path_tmp" ::: {params.samples}
        # Compress each sample
        ls $path_breast/*.frags.tsv | xargs -n 1 -P {threads} bash workflow/scripts/dts/breast/compress_sample.sh && \
        rm $path_breast/tmp.tsv.bgz
        """
        
        
rule prcannot_breast:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        img='workflow/envs/gretabench.sif',
        gid=rules.gen_gid_ensmbl.output.hg38,
    output:
        rna=temp(local('dts/hg38/breast/rna.h5ad')),
        ann=temp(local('dts/hg38/breast/annot.csv'))
    params:
        adata=config['dts']['breast']['url']['anndata'],
    shell:
        """
        path_breast=$(dirname {output.rna})
        wget --no-verbose '{params.adata}' -O "${{path_breast}}/tmp.h5ad"
        python workflow/scripts/dts/breast/annot.py \
        "${{path_breast}}/tmp.h5ad" \
        {input.gid} \
        {output.rna} \
        {output.ann}
        rm $path_breast/tmp.h5ad
        """


rule callpeaks_breast:
    threads: 8
    singularity: 'workflow/envs/gretabench.sif'
    input:
        frags=rules.download_fragments_breast.output.frag,
        annot=rules.prcannot_breast.output.ann,
    output: peaks=temp(local('dts/hg38/breast/peaks.h5ad'))
    resources:
        mem_mb=128000,
        runtime=720,
    shell:
        """
        python workflow/scripts/dts/callpeaks.py \
        -f {input.frags} \
        -a {input.annot} \
        -t $TMPDIR \
        -n {threads} \
        -o {output.peaks}
        """


rule annotate_breast:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        path_rna=rules.prcannot_breast.output.rna,
        path_peaks=rules.callpeaks_breast.output.peaks,
        path_ann=rules.prcannot_breast.output.ann,
    output: out='dts/hg38/breast/annotated.h5mu'
    shell:
        """
        python workflow/scripts/dts/breast/breast.py \
        -a {input.path_rna} \
        -b {input.path_peaks} \
        -c {input.path_ann} \
        -f {output}
        """
