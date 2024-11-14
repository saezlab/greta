rule download_pbmc10k:
    threads: 1
    singularity: 'workflow/envs/figr.sif'
    output:
        frags='dts/pbmc10k/smpl.frags.tsv.gz',
        tbis='dts/pbmc10k/smpl.frags.tsv.gz.tbi',
    params:
        matrix=config['dts']['pbmc10k']['url']['matrix'],
        atac_frags=config['dts']['pbmc10k']['url']['atac_frags'],
    shell:
        """
        wget --no-verbose '{params.atac_frags}' -O '{output.frags}'
        bash workflow/scripts/dts/format_frags.sh {output.frags}
        """


rule prcannot_pbmc10k:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    output: annot=temp(local('dts/pbmc10k/annot.csv')),
    shell:
        "python workflow/scripts/dts/pbmc10k/prc_annot.py -a {output.annot}"


rule callpeaks_pbmc10k:
    threads: 32
    singularity: 'workflow/envs/gretabench.sif'
    input:
        frags=rules.download_pbmc10k.output.frags,
        annot=rules.prcannot_pbmc10k.output.annot,
    output: peaks=temp(local('dts/pbmc10k/peaks.h5ad'))
    resources: mem_mb=64000
    shell:
        """
        python workflow/scripts/dts/callpeaks.py \
        -f {input.frags} \
        -a {input.annot} \
        -t '/tmp/pbcm10k/' \
        -n {threads} \
        -o {output.peaks}
        """


rule annotate_pbmc10k:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        annot=rules.prcannot_pbmc10k.output.annot,
        peaks=rules.callpeaks_pbmc10k.output.peaks,
        gid=rules.gen_gid_ensmbl.output,
    output: out='dts/pbmc10k/annotated.h5mu'
    resources: mem_mb=32000
    shell:
        """
        python workflow/scripts/dts/pbmc10k/pbmc10k.py \
        -b {input.annot} \
        -c {input.gid} \
        -e {input.peaks} \
        -f {output.out}
        """
