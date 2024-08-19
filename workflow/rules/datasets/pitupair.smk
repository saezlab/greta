rule download_pitupair:
    output:
        multi=temp(local('datasets/pitupair/multiome_original.h5')),
        frags=temp(local('datasets/pitupair/smpl.frags.tsv.gz'))
    params:
        multi=config['datasets']['pitupair']['url']['multi'],
        frags=config['datasets']['pitupair']['url']['frags'],
    shell:
        """
        wget '{params.frags}' -O '{output.frags}'
        wget '{params.multi}' -O '{output.multi}'
        """


rule prcannot_pitupair:
    input:
        'datasets/pitupair/multiome_original.h5'
    output:
        temp(local('datasets/pitupair/annot.csv'))
    shell:
        """
        python workflow/scripts/datasets/pitupair/prc_annot.py \
        -b {input} \
        -c {output}
        """


rule callpeaks_pitupair:
    input:
        frags='datasets/pitupair/smpl.frags.tsv.gz',
        annot='datasets/pitupair/annot.csv',
    output:
        tmp=temp(directory(local('datasets/pitupair/tmp_peaks'))),
        peaks=temp(local('datasets/pitupair/peaks.h5ad'))
    resources:
        mem_mb=64000,
    threads: 16
    shell:
        """
        python workflow/scripts/datasets/callpeaks.py \
        -f {input.frags} \
        -a {input.annot} \
        -t {output.tmp} \
        -o {output.peaks}
        """


rule annotate_pitupair:
    input:
        annot='datasets/pitupair/annot.csv',
        g='gdata/geneids',
        peaks='datasets/pitupair/peaks.h5ad',
        multi='datasets/pitupair/multiome_original.h5',
    output:
        'datasets/pitupair/annotated.h5mu'
    params:
        organism=config['datasets']['pitupair']['organism'],
    resources:
        mem_mb=32000,
    shell:
        """
        python workflow/scripts/datasets/pitupair/pitupair.py \
        -b {input.annot} \
        -c {input.g}
        -d {params.organism} \
        -e {input.peaks} \
        -f {output} \
        -g {input.multi}
        """
