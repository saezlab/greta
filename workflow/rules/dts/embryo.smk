localrules: download_embryo, annotate_embryo


rule download_embryo:
    threads: 3
    singularity: 'workflow/envs/figr.sif'
    input:
        img='workflow/envs/figr.sif',
        gid=rules.gen_gid_ensmbl.output.hg38,
    output:
        ann=temp(local('dts/hg38/embryo/annot.csv')),
        frag=temp(local(expand('dts/hg38/embryo/{sample}.frags.tsv.gz', sample=config['dts']['embryo']['samples']))),
        rna=temp(local('dts/hg38/embryo/rna.h5ad')),
    params:
        url_tar=config['dts']['embryo']['url']['tar'],
        url_rna=config['dts']['embryo']['url']['rna'],
        url_ann=config['dts']['embryo']['url']['ann'],
        samples=config['dts']['embryo']['samples'],
    shell:
        """
        path_embryo=$(dirname {output.ann})
        path_tmp=$path_embryo/tmp.tar
        path_rds=$path_embryo/rna.rds.gz
        path_h5=$path_embryo/rna.h5
        wget --no-verbose "{params.url_ann}" -O {output.ann}.gz
        gunzip {output.ann}
        wget --no-verbose "{params.url_tar}" -O $path_tmp
        wget --no-verbose "{params.url_rna}" -O $path_rds.gz
        gunzip $path_rds.gz
        tar -xvf $path_tmp -C $path_embryo
        rm $path_embryo/GSM6739292_*
        rm $path_embryo/GSM6739293_*
        rm $path_embryo/GSM6739294_*
        Rscript -e "
        library(rhdf5)
        library(Matrix)
        m <- readRDS('$path_rds')
        h5createFile('$path_h5')
        h5write(rownames(m), '$path_h5', 'var')
        h5write(colnames(m), '$path_h5', 'obs')
        h5write(m@i, '$path_h5', 'i')
        h5write(m@p, '$path_h5', 'p')
        h5write(m@x, '$path_h5', 'x')"
        python workflow/scripts/dts/embryo/prcrna.py \
        $path_h5 \
        {output.ann} \
        {input.gid}
        rm $path_rds
        rm $path_h5
        for f in $path_embryo/*_model_atac_fragments.tsv.gz; do
            [ -e "$f" ] || continue
            base="${{f%_model_atac_fragments.tsv.gz}}"
            idb="${{base#*_}}"
            new="${{idb}}.frags.tsv.gz"
            mv -- "$f" $path_embryo/$new
        done
        ls $path_embryo/*.frags.tsv.gz | xargs -n 1 -P {threads} bash workflow/scripts/dts/format_frags.sh
        rm $path_tmp
        """


rule callpeaks_embryo:
    threads: 8
    singularity: 'workflow/envs/gretabench.sif'
    input:
        frags=rules.download_embryo.output.frag,
        annot=rules.download_embryo.output.ann,
    output: peaks=temp(local('dts/hg38/embryo/peaks.h5ad'))
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


rule annotate_embryo:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        rna=rules.download_embryo.output.rna,
        peaks=rules.callpeaks_embryo.output.peaks,
        annot=rules.download_embryo.output.ann,
    output: out='dts/hg38/embryo/annotated.h5mu'
    shell:
        """
        python workflow/scripts/dts/embryo/embryo.py \
        -a {input.rna} \
        -b {input.peaks} \
        -c {input.annot} \
        -d {output.out}
        """
