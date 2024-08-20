rule download_pitunpair:
    output:
        gex=temp(local('datasets/pitunpair/smpl.filtered_feature_bc_matrix.h5')),
        peaks=temp(local('datasets/pitunpair/peaks.original.h5')),
        frags=temp(local('datasets/pitunpair/smpl.frags.tsv.gz')),
        fragIndex=temp(local('datasets/pitunpair/smpl.frags.tsv.gz.tbi'))
    params:
        gex=config['datasets']['pitunpair']['url']['rna_mtx'],
        peaks=config['datasets']['pitunpair']['url']['peaks'],
        frags=config['datasets']['pitunpair']['url']['atac_frags'],
        unzip='datasets/pitunpair/smpl.frags.tsv'
    shell:
        """
        wget '{params.gex}' -O '{output.gex}'
        wget '{params.peaks}' -O '{output.peaks}'
        wget '{params.frags}' -O '{output.frags}'
        gunzip -d {output.frags}
        bgzip {params.unzip}
        tabix -p bed {output.frags}
        """


rule coembedd_pitunpair:
    input:
        gex='datasets/pitunpair/smpl.filtered_feature_bc_matrix.h5',
        peaks='datasets/pitunpair/peaks.original.h5',
        frags='datasets/pitunpair/smpl.frags.tsv.gz'
    output:
        annot=temp(local('datasets/pitunpair/annot.csv')),
        exprMat=temp(local('datasets/pitunpair/exprMat.rds')),
        atacSE=temp(local('datasets/pitunpair/atac.se.rds')),
        cca=temp(local('datasets/pitunpair/cca.rds'))
    singularity: 
        'workflow/envs/figr.sif'
    shell:
        """
        Rscript workflow/scripts/datasets/pitunpair/coembedd.R \
        {input.gex} \
        {input.peaks} \
        {input.frags} \
        {output.annot} \
        {output.exprMat} \
        {output.atacSE} \
        {output.cca} 
        """


rule paircells_pitunpair:
    input:
        exprMat='datasets/pitunpair/exprMat.rds',
        atacSE='datasets/pitunpair/atac.se.rds',
        cca='datasets/pitunpair/cca.rds'
    output:
        barMap=temp(local('datasets/pitunpair/barMap.csv'))
    singularity:
        'workflow/envs/figr.sif'
    shell:
        """
        Rscript workflow/scripts/datasets/pitunpair/pairCells.R \
        {input.exprMat} \
        {input.atacSE} \
        {input.cca} \
        {output.barMap}
        """


rule callpeaks_pitunpair:
    threads: 32
    input:
        frags='datasets/pitunpair/smpl.frags.tsv.gz',
        annot='datasets/pitunpair/annot.csv',
    output:
        tmp=temp(directory(local('datasets/pitunpair/tmp_peaks'))),
        peaks=temp(local('datasets/pitunpair/peaks.h5ad'))
    resources:
        mem_mb=32000,
    shell:
        """
        python workflow/scripts/datasets/callpeaks.py \
        -f {input.frags} \
        -a {input.annot} \
        -t {output.tmp} \
        -o {output.peaks}
        """


rule annotate_pitunpair:
    input:
        annot='datasets/pitunpair/annot.csv',
        g='gdata/geneids',
        peaks='datasets/pitunpair/peaks.h5ad',
        gex='datasets/pitunpair/smpl.filtered_feature_bc_matrix.h5',
        barmap='datasets/pitunpair/barMap.csv'
    output:
        tmp=temp(directory(local('datasets/pitunpair/tmp_annot'))),
        out='datasets/pitunpair/annotated.h5mu'
    params:
        organism=config['datasets']['pitunpair']['organism'],
    resources:
        mem_mb=32000,
    shell:
        """
        python workflow/scripts/datasets/pitunpair/pitunpair.py \
        -a {output.tmp} \
        -b {input.annot} \
        -c {input.g} \
        -d {params.organism} \
        -e {input.peaks} \
        -f {output.out} \
        -g {input.gex} \
        -i {input.barmap}
        """
