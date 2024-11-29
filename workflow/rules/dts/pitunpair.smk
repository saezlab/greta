rule download_pitunpair:
    threads: 1
    singularity: 'workflow/envs/figr.sif'
    output:
        gex=temp(local('dts/pitunpair/smpl.filtered_feature_bc_matrix.h5')),
        peaks=temp(local('dts/pitunpair/peaks.original.h5')),
        frags='dts/pitunpair/smpl.frags.tsv.gz',
        tbis='dts/pitunpair/smpl.frags.tsv.gz.tbi',
        annot=temp(local('dts/pitunpair/annot.csv')),
    params:
        gex=config['dts']['pitunpair']['url']['rna_mtx'],
        peaks=config['dts']['pitunpair']['url']['peaks'],
        frags=config['dts']['pitunpair']['url']['atac_frags'],
        annot=config['dts']['pitunpair']['url']['annot']
    shell:
        """
        wget --no-verbose '{params.gex}' -O '{output.gex}'
        wget --no-verbose '{params.peaks}' -O '{output.peaks}'
        wget --no-verbose '{params.frags}' -O '{output.frags}'
        bash workflow/scripts/dts/format_frags.sh {output.frags}
        wget --no-verbose '{params.annot}' -O '{output.annot}'
        """


rule coembedd_pitunpair:
    threads: 32
    singularity: 'workflow/envs/figr.sif'
    input:
        gex=rules.download_pitunpair.output.gex,
        annot=rules.download_pitunpair.output.annot,
        peaks=rules.download_pitunpair.output.peaks,
        frags=rules.download_pitunpair.output.frags,
        tbis=rules.download_pitunpair.output.tbis,
    output: cca=temp(local('dts/pitunpair/cca.rds'))
    shell:
        """
        Rscript workflow/scripts/dts/pitunpair/coembedd.R \
        {input.gex} \
        {input.annot} \
        {input.peaks} \
        {input.frags} \
        {output.cca}
        """


rule paircells_pitunpair:
    threads: 1
    singularity: 'workflow/envs/figr.sif'
    input:
        cca=rules.coembedd_pitunpair.output.cca,
        annot=rules.download_pitunpair.output.annot,
    output: barmap=temp(local('dts/pitunpair/barmap.csv'))
    shell:
        """
        Rscript workflow/scripts/dts/pitunpair/paircells.R \
        {input.cca} \
        {input.annot} \
        {output.barmap}
        """


rule callpeaks_pitunpair:
    threads: 32
    singularity: 'workflow/envs/gretabench.sif'
    input:
        frags=rules.download_pitunpair.output.frags,
        annot=rules.paircells_pitunpair.output.barmap,
    output: peaks=temp(local('dts/pitunpair/peaks.h5ad')),
    resources: mem_mb=32000
    shell:
        """
        python workflow/scripts/dts/callpeaks.py \
        -f {input.frags} \
        -a {input.annot} \
        -t '/tmp/pitunpair/' \
        -n {threads} \
        -o {output.peaks}
        """


rule annotate_pitunpair:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        peaks=rules.callpeaks_pitunpair.output.peaks,
        gex=rules.download_pitunpair.output.gex,
        barmap=rules.paircells_pitunpair.output.barmap,
        gid=rules.gen_gid_ensmbl.output,
    output:
        out='dts/pitunpair/annotated.h5mu'
    resources: mem_mb=32000
    shell:
        """
        python workflow/scripts/dts/pitunpair/pitunpair.py \
        -c {input.gid} \
        -e {input.peaks} \
        -f {output.out} \
        -g {input.gex} \
        -i {input.barmap}
        """

