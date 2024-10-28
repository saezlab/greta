rule download_pitunpair:
    threads: 1
    singularity: 'workflow/envs/figr.sif'
    output:
        gex=temp(local('datasets/pitunpair/smpl.filtered_feature_bc_matrix.h5')),
        peaks=temp(local('datasets/pitunpair/peaks.original.h5')),
        frags='datasets/pitunpair/smpl.frags.tsv.gz',
        tbis='datasets/pitunpair/smpl.frags.tsv.gz.tbi',
        annot=temp(local('datasets/pitunpair/annot.csv')),
    params:
        gex=config['datasets']['pitunpair']['url']['rna_mtx'],
        peaks=config['datasets']['pitunpair']['url']['peaks'],
        frags=config['datasets']['pitunpair']['url']['atac_frags'],
        annot=config['datasets']['pitunpair']['url']['annot']
    shell:
        """
        wget --no-verbose '{params.gex}' -O '{output.gex}'
        wget --no-verbose '{params.peaks}' -O '{output.peaks}'
        wget --no-verbose '{params.frags}' -O '{output.frags}'
        bash workflow/scripts/datasets/format_frags.sh {output.frags}
        wget --no-verbose '{params.annot}' -O '{output.annot}'
        """


rule coembedd_pitunpair:
    threads: 32
    input:
        gex=rules.download_pitunpair.output.gex,
        annot=rules.download_pitunpair.output.annot,
        peaks=rules.download_pitunpair.output.peaks,
        frags=rules.download_pitunpair.output.frags,
        tbis=rules.download_pitunpair.output.tbis,
    output:
        cca=temp(local('datasets/pitunpair/cca.rds'))
    singularity: 
        'workflow/envs/figr.sif'
    shell:
        """
        Rscript workflow/scripts/datasets/pitunpair/coembedd.R \
        {input.gex} \
        {input.annot} \
        {input.peaks} \
        {input.frags} \
        {output.cca}
        """


rule paircells_pitunpair:
    threads: 1
    input:
        cca=rules.coembedd_pitunpair.output.cca,
        annot=rules.download_pitunpair.output.annot,
    output:
        barmap=temp(local('datasets/pitunpair/barmap.csv'))
    singularity:
        'workflow/envs/figr.sif'
    shell:
        """
        Rscript workflow/scripts/datasets/pitunpair/paircells.R \
        {input.cca} \
        {input.annot} \
        {output.barmap}
        """


rule callpeaks_pitunpair:
    threads: 32
    input:
        frags=rules.download_pitunpair.output.frags,
        annot=rules.paircells_pitunpair.output.barmap,
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        tmp=temp(directory(local('datasets/pitunpair/tmp_peaks'))),
        peaks=temp(local('datasets/pitunpair/peaks.h5ad')),
    resources:
        mem_mb=32000,
    shell:
        """
        python workflow/scripts/datasets/callpeaks.py \
        -f {input.frags} \
        -a {input.annot} \
        -t {output.tmp} \
        -n {threads} \
        -o {output.peaks}
        """


rule annotate_pitunpair:
    threads: 1
    input:
        peaks=rules.callpeaks_pitunpair.output.peaks,
        gex=rules.download_pitunpair.output.gex,
        barmap=rules.paircells_pitunpair.output.barmap,
        g=rules.download_geneids.output.dr,
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        out='datasets/pitunpair/annotated.h5mu'
    params:
        organism=config['datasets']['pitunpair']['organism'],
    resources:
        mem_mb=32000,
    shell:
        """
        python workflow/scripts/datasets/pitunpair/pitunpair.py \
        -c {input.g} \
        -d {params.organism} \
        -e {input.peaks} \
        -f {output.out} \
        -g {input.gex} \
        -i {input.barmap}
        """

