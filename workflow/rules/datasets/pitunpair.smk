rule download_pitunpair:
    threads: 1
    output:
        gex=temp(local('datasets/pitunpair/smpl.filtered_feature_bc_matrix.h5')),
        peaks=temp(local('datasets/pitunpair/peaks.original.h5')),
        frags='datasets/pitunpair/smpl.frags.tsv.gz',
        celltypes=temp(local('datasets/pitunpair/celltypes.csv')),
    params:
        gex=config['datasets']['pitunpair']['url']['rna_mtx'],
        peaks=config['datasets']['pitunpair']['url']['peaks'],
        frags=config['datasets']['pitunpair']['url']['atac_frags'],
        celltypes=config['datasets']['pitunpair']['url']['celltypes']
    shell:
        """
        wget --no-verbose '{params.gex}' -O '{output.gex}'
        wget --no-verbose '{params.peaks}' -O '{output.peaks}'
        wget --no-verbose '{params.frags}' -O '{output.frags}'
        bash workflow/scripts/datasets/format_frags.sh {output.frags}
        wget --no-verbose '{params.celltypes}' -O '{output.celltypes}'
        """


rule index_frags:
    threads: 1
    input:
        frags=rules.download_pitunpair.output.frags
    output:
        index=temp(local('datasets/pitunpair/smpl.frags.tsv.gz.tbi'))
    params:
        unzip=rules.download_pitunpair.output.frags.replace('.gz', '')
    singularity:
        'workflow/envs/figr.sif'
    shell:
        """
        gunzip -d {input.frags}
        bgzip {params.unzip}
        tabix -p bed {input.frags}
        """


rule coembedd_pitunpair:
    threads: 32
    input:
        gex=rules.download_pitunpair.output.gex,
        celltypes=rules.download_pitunpair.output.celltypes,
        peaks=rules.download_pitunpair.output.peaks,
        frags=rules.download_pitunpair.output.frags,
        index=rules.index_frags.output.index
    output:
        cca=temp(local('datasets/pitunpair/cca.rds'))
    singularity: 
        'workflow/envs/figr.sif'
    shell:
        """
        Rscript workflow/scripts/datasets/pitunpair/coembedd.R \
        {input.gex} \
        {input.celltypes} \
        {input.peaks} \
        {input.frags} \
        {output.cca} \
        {threads}
        """


rule paircells_pitunpair:
    threads: 1
    input:
        cca=rules.coembedd_pitunpair.output.cca,
        celltypes=rules.download_pitunpair.output.celltypes,
    output:
        barmap=temp(local('datasets/pitunpair/barmap.csv'))
    singularity:
        'workflow/envs/figr.sif'
    shell:
        """
        Rscript workflow/scripts/datasets/pitunpair/paircells.R \
        {input.cca} \
        {input.celltypes} \
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

