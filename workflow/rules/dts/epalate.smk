localrules: download_epalate, annotate_epalate, download_epalate_annot


rule download_epalate_frag:
    threads: 2
    singularity: 'workflow/envs/gretabench.sif'
    input:
        img='workflow/envs/gretabench.sif',
    output:
        frag=local(expand(
            'dts/mm10/epalate/.tmp/{sample}.frags.tsv.gz',
            sample=config['dts']['epalate']['samples'],
        )),
        tbi=local(expand(
            'dts/mm10/epalate/.tmp/{sample}.frags.tsv.gz.tbi',
            sample=config['dts']['epalate']['samples'],
        )),
    params:
        url_frag_E1401=config['dts']['epalate']['url']['frag_E1401'],
        url_frag_E1251=config['dts']['epalate']['url']['frag_E1251'],
    shell:
        """
        mkdir -p dts/mm10/epalate/.tmp

        wget --no-verbose --continue --tries=5 --waitretry=2 --content-disposition \
             '{params.url_frag_E1401}' \
             -O dts/mm10/epalate/.tmp/E1401.frags.tsv.gz

        wget --no-verbose --continue --tries=5 --waitretry=2 --content-disposition \
             '{params.url_frag_E1251}' \
             -O dts/mm10/epalate/.tmp/E1251.frags.tsv.gz

        ls dts/mm10/epalate/.tmp/*.frags.tsv.gz | xargs -n 1 -P {threads} bash workflow/scripts/dts/format_frags.sh
        """


rule download_epalate_annot:
    threads: 2
    singularity: 'workflow/envs/figr.sif'
    input:
        img='workflow/envs/figr.sif',
    output:
        annot=temp(local('dts/mm10/epalate/.tmp/annot.csv')),
    params:
        url_gdrive=config['dts']['epalate']['url']['gdrive'],
    shell:
        """
        mkdir -p dts/mm10/epalate/.tmp
        rds_tmp=dts/mm10/epalate/.tmp/processed.rds

        wget --no-verbose --continue --tries=5 --waitretry=2 --content-disposition \
             '{params.url_gdrive}' \
             -O "$rds_tmp"

        Rscript workflow/scripts/dts/epalate/rds_to_annot.R \
            "$rds_tmp" \
            "{output.annot}"
        """


rule download_epalate_rna:
    threads: 2
    singularity: 'workflow/envs/gretabench.sif'
    input:
        img='workflow/envs/gretabench.sif',
        ann='dts/mm10/epalate/.tmp/annot.csv',
    output:
        rna_E1401=temp(local('dts/mm10/epalate/.tmp/E1401_rna.h5ad')),
        rna_E1251=temp(local('dts/mm10/epalate/.tmp/E1251_rna.h5ad')),
    params:
        url_h5_E1401=config['dts']['epalate']['url']['h5_E1401'],
        url_h5_E1251=config['dts']['epalate']['url']['h5_E1251'],
    shell:
        """
        mkdir -p dts/mm10/epalate/.tmp

        h5_tmp_E1401=dts/mm10/epalate/.tmp/E1401_filtered_feature_bc_matrix.h5
        wget --no-verbose --continue --tries=5 --waitretry=2 --content-disposition \
             '{params.url_h5_E1401}' \
             -O "$h5_tmp_E1401"
        python workflow/scripts/dts/epalate/h5_to_h5ad.py \
            -a "$h5_tmp_E1401" \
            -b {input.ann} \
            -c E1401 \
            -d {output.rna_E1401}

        h5_tmp_E1251=dts/mm10/epalate/.tmp/E1251_filtered_feature_bc_matrix.h5
        wget --no-verbose --continue --tries=5 --waitretry=2 --content-disposition \
             '{params.url_h5_E1251}' \
             -O "$h5_tmp_E1251"
        python workflow/scripts/dts/epalate/h5_to_h5ad.py \
            -a "$h5_tmp_E1251" \
            -b {input.ann} \
            -c E1251 \
            -d {output.rna_E1251}
        """


rule download_epalate:
    threads: 2
    singularity: 'workflow/envs/gretabench.sif'
    input:
        img='workflow/envs/gretabench.sif',
        frag=rules.download_epalate_frag.output.frag,
        frag_tbi=rules.download_epalate_frag.output.tbi,
        rna_E1401=rules.download_epalate_rna.output.rna_E1401,
        rna_E1251=rules.download_epalate_rna.output.rna_E1251,
        annot=rules.download_epalate_annot.output.annot,
    output:
        frag=local(expand(
            'dts/mm10/epalate/{sample}.frags.tsv.gz',
            sample=config['dts']['epalate']['samples'],
        )),
        tbi=local(expand(
            'dts/mm10/epalate/{sample}.frags.tsv.gz.tbi',
            sample=config['dts']['epalate']['samples'],
        )),
        rna=local('dts/mm10/epalate/rna.h5ad'),
        annot=local('dts/mm10/epalate/annot.csv'),
    shell:
        """
        python workflow/scripts/dts/epalate/filter_merge.py \
        -a {input.frag} \
        -b {input.rna_E1401} {input.rna_E1251} \
        -c {input.annot} \
        -d {output.frag} \
        -e {output.rna} \
        -f {output.annot}
        rm -rf dts/mm10/epalate/.tmp
        """


rule callpeaks_epalate:
    threads: 8
    singularity: 'workflow/envs/gretabench.sif'
    input:
        img='workflow/envs/gretabench.sif',
        frag=rules.download_epalate.output.frag,
        annot=rules.download_epalate.output.annot,
    output: peaks=temp(local('dts/mm10/epalate/peaks.h5ad'))
    resources:
        mem_mb=128000,
        runtime=2160,
    shell:
        """
        python workflow/scripts/dts/callpeaks_mm10.py \
        -f {input.frag} \
        -a {input.annot} \
        -t $TMPDIR \
        -n {threads} \
        -o {output.peaks}
        """


rule annotate_epalate:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        img='workflow/envs/gretabench.sif',
        rna=rules.download_epalate.output.rna,
        peaks=rules.callpeaks_epalate.output.peaks,
        annot=rules.download_epalate.output.annot,
    output: out='dts/mm10/epalate/annotated.h5mu'
    shell:
        """
        python workflow/scripts/dts/epalate/epalate.py \
        -a {input.rna} \
        -b {input.peaks} \
        -c {input.annot} \
        -d {output.out}
        """
