rule download_pitupair:
    output:
        multi=temp(local('datasets/pitupaired/multiome_original.h5')),
        frags=temp(local('datasets/pitupaired/smpl.frags.tsv.gz'))
    params:
        multi=config['datasets']['pitupaired']['url']['multi'],
        frags=config['datasets']['pitupaired']['url']['frags'],
    shell:
        """
        wget '{params.frags}' -O '{output.frags}'
        wget '{params.multi}' -O '{output.multi}'
        """


rule prcannot_pitupair:
    input:
        multi='datasets/pitupaired/multiome_original.h5'
    output:
        tmp=temp(directory(local('datasets/pitupaired/tmp'))),
        annot=temp(local('datasets/pitupaired/annot.csv'))
    shell:
        """
        python workflow/scripts/datasets/pitupaired/prc_annot.py \
        -a {output.tmp} \
        -b {input.multi} \
        -c {output.annot}
        """


rule callpeaks_pitupair:
    input:
        frags='datasets/pitupaired/smpl.frags.tsv.gz',
        annot='datasets/pitupaired/annot.csv',
    output:
        tmp=temp(directory(local('datasets/pitupaired/tmp_peaks'))),
        peaks=temp(local('datasets/pitupaired/peaks.h5ad'))
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
        annot='datasets/pitupaired/annot.csv',
        g='gdata/geneids',
        peaks='datasets/pitupaired/peaks.h5ad',
        multi='datasets/pitupaired/multiome_original.h5',
    output:
        tmp=temp(directory(local('datasets/pitupaired/tmp_annot'))),
        out=temp('datasets/pitupaired/annotated.h5mu')
    params:
        organism=config['datasets']['pitupaired']['organism'],
    resources:
        mem_mb=32000,
    shell:
        """
        python workflow/scripts/datasets/pitupaired/pitupaired.py \
        -a {output.tmp} \
        -b {input.annot} \
        -c {input.g}
        -d {params.organism} \
        -e {input.peaks} \
        -f {output.out} \
        -g {input.multi}
        """
