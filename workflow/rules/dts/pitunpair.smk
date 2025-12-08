rule download_pitunpair:
    threads: 1
    singularity: 'workflow/envs/figr.sif'
    input: 'workflow/envs/figr.sif'
    output:
        gex=temp(local('dts/hg38/pitunpair/smpl.filtered_feature_bc_matrix.h5')),
        peaks=temp(local('dts/hg38/pitunpair/peaks.original.h5')),
        frags='dts/hg38/pitunpair/smpl.frags.tsv.gz',
        tbis='dts/hg38/pitunpair/smpl.frags.tsv.gz.tbi',
        annot=temp(local('dts/hg38/pitunpair/annot.csv')),
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
    threads: 4
    conda: '../../envs/glue.yaml'
    input:
        gex=rules.download_pitunpair.output.gex,
        acc=rules.download_pitunpair.output.peaks,
        gid=rules.gen_gid_ensmbl.output.hg38,
        ann=rules.download_pitunpair.output.annot,
    output: temp(local('dts/hg38/pitunpair/X_glue.csv'))
    resources:
        partition='gpu-single',
        slurm='gres=gpu:1'
    shell:
        """
        python workflow/scripts/dts/coembed.py \
        {input.gex} \
        {input.acc} \
        {input.gid} \
        {input.ann} \
        {output}
        """


rule paircells_pitunpair:
    threads: 1
    singularity: 'workflow/envs/figr.sif'
    input:
        emb=rules.coembedd_pitunpair.output,
        ann=rules.download_pitunpair.output.annot,
    output: temp(local('dts/hg38/pitunpair/barmap.csv'))
    shell:
        """
        Rscript workflow/scripts/dts/pitunpair/paircells.R \
        {input.emb} \
        {input.ann} \
        {output}
        """


rule callpeaks_pitunpair:
    threads: 32
    singularity: 'workflow/envs/gretabench.sif'
    input:
        img='workflow/envs/gretabench.sif',
        frags=rules.download_pitunpair.output.frags,
        annot=rules.paircells_pitunpair.output,
    output: peaks=temp(local('dts/hg38/pitunpair/peaks.h5ad')),
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
        barmap=rules.paircells_pitunpair.output,
        gid=rules.gen_gid_ensmbl.output.hg38,
    output:
        out='dts/hg38/pitunpair/annotated.h5mu'
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

