localrules: download_geneids, download_pbmc10k


rule download_geneids:
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        hg='gdata/geneids/hg38.csv',
        mm='gdata/geneids/mm10.csv',
        dr=directory('gdata/geneids')
    shell:
        """
        Rscript workflow/scripts/datasets/download_geneids.R \
        {output.hg} \
        {output.mm}
        """


rule download_pbmc10k:
    output:
        atac_frags='datasets/pbmc10k/smpl.frags.tsv.gz',
    params:
        matrix=config['datasets']['pbmc10k']['url']['matrix'],
        atac_frags=config['datasets']['pbmc10k']['url']['atac_frags'],
    shell:
        """
        wget --no-verbose '{params.atac_frags}' -O '{output.atac_frags}'
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
    threads: 32
    input:
        frags=rules.download_pbmc10k.output.atac_frags,
        annot=rules.prcannot_pbmc10k.output.annot,
    output:
        tmp=temp(directory(local('datasets/pbmc10k/tmp_peaks'))),
        peaks=temp(local('datasets/pbmc10k/peaks.h5ad'))
    resources:
        mem_mb=64000,
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
        annot=rules.prcannot_pbmc10k.output.annot,
        g=rules.download_geneids.output.dr,
        peaks=rules.callpeaks_pbmc10k.output.peaks,
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
