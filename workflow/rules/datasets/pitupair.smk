rule download_pitupair:
    threads: 1
    singularity: 'workflow/envs/figr.sif'
    output:
        gex=temp(local('datasets/pitupair/multiome_original.h5')),
        frags='datasets/pitupair/smpl.frags.tsv.gz',
        tbis='datasets/pitupair/smpl.frags.tsv.gz.tbi',
        annot=temp(local('datasets/pitupair/annot.csv'))
    params:
        gex=config['datasets']['pitupair']['url']['gex'],
        frags=config['datasets']['pitupair']['url']['frags'],
        annot=config['datasets']['pitupair']['url']['annot']
    
    shell:
        """
        wget --no-verbose '{params.frags}' -O '{output.frags}'
        bash workflow/scripts/datasets/format_frags.sh {output.frags}
        wget --no-verbose '{params.gex}' -O '{output.gex}'
        wget --no-verbose '{params.annot}' -O '{output.annot}'
        awk 'BEGIN {{FS=OFS=","}} NR==1 {{print $0; next}} {{gsub(/-[0-9]+$/, "", $1); print $3"_"$1,$2,$3}}' {output.annot} > {output.annot}.tmp
        mv {output.annot}.tmp {output.annot}
        """


rule callpeaks_pitupair:
    threads: 32
    input:
        frags=rules.download_pitupair.output.frags,
        annot=rules.download_pitupair.output.annot,
    output:
        tmp=temp(directory(local('datasets/pitupair/tmp_peaks'))),
        peaks=temp(local('datasets/pitupair/peaks.h5ad'))
    singularity:
        'workflow/envs/gretabench.sif'
    resources:
        mem_mb=64000,
    shell:
        """
        python workflow/scripts/datasets/callpeaks.py \
        -f {input.frags} \
        -a {input.annot} \
        -t {output.tmp} \
        -n {threads} \
        -o {output.peaks}
        """


rule annotate_pitupair:
    threads: 1
    input:
        annot=rules.download_pitupair.output.annot,
        peaks=rules.callpeaks_pitupair.output.peaks,
        gex=rules.download_pitupair.output.gex,
        g=rules.download_geneids.output.dr,
    output:
        out='datasets/pitupair/annotated.h5mu'
    singularity:
        'workflow/envs/gretabench.sif'
    params:
        organism=config['datasets']['pitupair']['organism'],
    resources:
        mem_mb=32000,
    shell:
        """
        python workflow/scripts/datasets/pitupair/pitupair.py \
        -b {input.annot} \
        -c {input.g} \
        -d {params.organism} \
        -e {input.peaks} \
        -f {output} \
        -g {input.gex}
        """
