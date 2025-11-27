localrules: download_mouse, annotate_mouse, download_mouse_annot


rule download_mouse_frag:
    threads: 2
    singularity: 'workflow/envs/gretabench.sif'
    input:
        img='workflow/envs/gretabench.sif',
    output:
        frag=local(expand(
            'dts/mm10/mouse/.tmp/{sample}.frags.tsv.gz',
            sample=config['dts']['mouse']['samples'],
        )),
        tbi=local(expand(
            'dts/mm10/mouse/.tmp/{sample}.frags.tsv.gz.tbi',
            sample=config['dts']['mouse']['samples'],
        )),
    params:
        url_frag_E1401=config['dts']['mouse']['url']['frag_E1401'],
        url_frag_E1251=config['dts']['mouse']['url']['frag_E1251'],
    shell:
        """
        mkdir -p dts/mm10/mouse/.tmp

        wget --no-verbose --continue --tries=5 --waitretry=2 --content-disposition \
             '{params.url_frag_E1401}' \
             -O dts/mm10/mouse/.tmp/E1401.frags.tsv.gz

        wget --no-verbose --continue --tries=5 --waitretry=2 --content-disposition \
             '{params.url_frag_E1251}' \
             -O dts/mm10/mouse/.tmp/E1251.frags.tsv.gz

        ls dts/mm10/mouse/.tmp/*.frags.tsv.gz | xargs -n 1 -P {threads} bash workflow/scripts/dts/format_frags.sh
        """


rule download_mouse_annot:
    threads: 2
    singularity: 'workflow/envs/figr.sif'
    input:
        img='workflow/envs/figr.sif',
    output:
        annot=temp(local('dts/mm10/mouse/.tmp/annot.csv')),
    params:
        url_gdrive=config['dts']['mouse']['url']['gdrive'],
    shell:
        """
        mkdir -p dts/mm10/mouse/.tmp
        rds_tmp=dts/mm10/mouse/.tmp/processed.rds

        wget --no-verbose --continue --tries=5 --waitretry=2 --content-disposition \
             '{params.url_gdrive}' \
             -O "$rds_tmp"

        Rscript workflow/scripts/dts/mouse/rds_to_annot.R \
            "$rds_tmp" \
            "{output.annot}"
        """


rule download_mouse_rna:
    threads: 2
    singularity: 'workflow/envs/gretabench.sif'
    input:
        img='workflow/envs/gretabench.sif',
        ann='dts/mm10/mouse/.tmp/annot.csv',
    output:
        rna_E1401=temp(local('dts/mm10/mouse/.tmp/E1401_rna.h5ad')),
        rna_E1251=temp(local('dts/mm10/mouse/.tmp/E1251_rna.h5ad')),
    params:
        url_h5_E1401=config['dts']['mouse']['url']['h5_E1401'],
        url_h5_E1251=config['dts']['mouse']['url']['h5_E1251'],
    shell:
        """
        mkdir -p dts/mm10/mouse/.tmp

        h5_tmp_E1401=dts/mm10/mouse/.tmp/E1401_filtered_feature_bc_matrix.h5
        wget --no-verbose --continue --tries=5 --waitretry=2 --content-disposition \
             '{params.url_h5_E1401}' \
             -O "$h5_tmp_E1401"
        python workflow/scripts/dts/mouse/h5_to_h5ad.py \
            -a "$h5_tmp_E1401" \
            -b {input.ann} \
            -c E1401 \
            -d {output.rna_E1401}

        h5_tmp_E1251=dts/mm10/mouse/.tmp/E1251_filtered_feature_bc_matrix.h5
        wget --no-verbose --continue --tries=5 --waitretry=2 --content-disposition \
             '{params.url_h5_E1251}' \
             -O "$h5_tmp_E1251"
        python workflow/scripts/dts/mouse/h5_to_h5ad.py \
            -a "$h5_tmp_E1251" \
            -b {input.ann} \
            -c E1251 \
            -d {output.rna_E1251}
        """


rule download_mouse:
    threads: 2
    singularity: 'workflow/envs/gretabench.sif'
    input:
        img='workflow/envs/gretabench.sif',
        frag=rules.download_mouse_frag.output.frag,
        frag_tbi=rules.download_mouse_frag.output.tbi,
        rna_E1401=rules.download_mouse_rna.output.rna_E1401,
        rna_E1251=rules.download_mouse_rna.output.rna_E1251,
        annot=rules.download_mouse_annot.output.annot,
    output:
        frag=local(expand(
            'dts/mm10/mouse/{sample}.frags.tsv.gz',
            sample=config['dts']['mouse']['samples'],
        )),
        tbi=local(expand(
            'dts/mm10/mouse/{sample}.frags.tsv.gz.tbi',
            sample=config['dts']['mouse']['samples'],
        )),
        rna=local('dts/mm10/mouse/rna.h5ad'),
        annot=local('dts/mm10/mouse/annot.csv'),
    shell:
        """
        python workflow/scripts/dts/mouse/filter_merge.py \
        -a {input.frag} \
        -b {input.rna_E1401} {input.rna_E1251} \
        -c {input.annot} \
        -d {output.frag} \
        -e {output.rna} \
        -f {output.annot}
        rm -rf dts/mm10/mouse/.tmp
        """


rule callpeaks_mouse:
    threads: 8
    singularity: 'workflow/envs/gretabench.sif'
    input:
        img='workflow/envs/gretabench.sif',
        frag=rules.download_mouse.output.frag,
        annot=rules.download_mouse.output.annot,
    output: peaks=temp(local('dts/mm10/mouse/peaks.h5ad'))
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


rule annotate_mouse:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        img='workflow/envs/gretabench.sif',
        rna=rules.download_mouse.output.rna,
        peaks=rules.callpeaks_mouse.output.peaks,
        annot=rules.download_mouse.output.annot,
    output: out='dts/mm10/mouse/annotated.h5mu'
    shell:
        """
        python workflow/scripts/dts/mouse/mouse.py \
        -a {input.rna} \
        -b {input.peaks} \
        -c {input.annot} \
        -d {output.out}
        """
