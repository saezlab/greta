rule download_pbmc10k:
    output:
        atac_frags='datasets/pbmc10k/smpl.frags.tsv.gz',
    params:
        matrix=config['datasets']['pbmc10k']['url']['matrix'],
        atac_frags=config['datasets']['pbmc10k']['url']['atac_frags'],
    shell:
        """
        wget '{params.atac_frags}' -O '{output.atac_frags}'
        """


rule prcannot_pbmc10k:
    output:
        tmp=temp(directory(local('datasets/pbmc10k/tmp'))),
        annot=temp(local('datasets/pbmc10k/annot.csv')),
    shell:
        """
        python workflow/scripts/datasets/pbmc10k/prc_annot.py \
        -t {output.tmp} \
        -a {output.annot}
        """


rule callpeaks_pbmc10k:
    input:
        frags='datasets/pbmc10k/smpl.frags.tsv.gz',
        annot='datasets/pbmc10k/annot.csv',
    output:
        tmp=temp(directory(local('datasets/pbmc10k/tmp_peaks'))),
        peaks=temp(local('datasets/pbmc10k/peaks.h5ad'))
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


rule annotate_pbmc10k:
    input:
        annot='datasets/pbmc10k/annot.csv',
        g='gdata/geneids',
        peaks='datasets/pbmc10k/peaks.h5ad',
    output:
        tmp=temp(directory(local('datasets/pbmc10k/tmp_annot'))),
        out='datasets/pbmc10k/annotated.h5mu'
    params:
        organism=config['datasets']['pbmc10k']['organism'],
    resources:
        mem_mb=32000,
    shell:
        """
        python workflow/scripts/datasets/pbmc10k/pbmc10k.py \
        -a {output.tmp} \
        -b {input.annot} \
        -c {input.g} \
        -d {params.organism} \
        -e {input.peaks} \
        -f {output.out}
        """
