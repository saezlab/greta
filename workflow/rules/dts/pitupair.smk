rule download_pitupair:
    threads: 1
    singularity: 'workflow/envs/figr.sif'
    output:
        gex=temp(local('dts/pitupair/multiome_original.h5')),
        frags='dts/pitupair/smpl.frags.tsv.gz',
        tbis='dts/pitupair/smpl.frags.tsv.gz.tbi',
        annot=temp(local('dts/pitupair/annot.csv'))
    params:
        gex=config['dts']['pitupair']['url']['gex'],
        frags=config['dts']['pitupair']['url']['frags'],
        annot=config['dts']['pitupair']['url']['annot']
    shell:
        """
        wget --no-verbose '{params.frags}' -O '{output.frags}' && \
        bash workflow/scripts/dts/format_frags.sh {output.frags} && \
        wget --no-verbose '{params.gex}' -O '{output.gex}' && \
        wget --no-verbose '{params.annot}' -O '{output.annot}' && \
        awk 'BEGIN {{FS=OFS=","}} NR==1 {{print $0; next}} {{gsub(/-[0-9]+$/, "", $1); print $3"_"$1,$2,$3}}' {output.annot} > {output.annot}.tmp && \
        mv {output.annot}.tmp {output.annot}
        """


rule callpeaks_pitupair:
    threads: 32
    singularity: 'workflow/envs/gretabench.sif'
    input:
        frags=rules.download_pitupair.output.frags,
        annot=rules.download_pitupair.output.annot,
    output: peaks=temp(local('dts/pitupair/peaks.h5ad'))
    resources: mem_mb=64000
    shell:
        """
        python workflow/scripts/dts/callpeaks.py \
        -f {input.frags} \
        -a {input.annot} \
        -t '/tmp/pitupair/' \
        -n {threads} \
        -o {output.peaks}
        """


rule annotate_pitupair:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        annot=rules.download_pitupair.output.annot,
        peaks=rules.callpeaks_pitupair.output.peaks,
        gex=rules.download_pitupair.output.gex,
        gid=rules.gen_gid_ensmbl.output,
    output: out='dts/pitupair/annotated.h5mu'
    resources: mem_mb=32000
    shell:
        """
        python workflow/scripts/dts/pitupair/pitupair.py \
        -b {input.annot} \
        -c {input.gid} \
        -e {input.peaks} \
        -f {output} \
        -g {input.gex}
        """
