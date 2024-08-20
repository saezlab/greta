rule download_pitunpaired:
    output:
        gex='datasets/pitunpaired/smpl.filtered_feature_bc_matrix.h5',
        peaks='datasets/pitunpaired/peaks.original.h5',
        frags='datasets/pitunpaired/smpl.frags.tsv.gz',
        fragIndex='datasets/pitunpaired/smpl.frags.tsv.gz.tbi'
    params:
        gex=config['datasets']['pitunpaired']['url']['rna_mtx'],
        peaks=config['datasets']['pitunpaired']['url']['peaks'],
        frags=config['datasets']['pitunpaired']['url']['atac_frags'],
        unzip='datasets/pitunpaired/smpl.frags.tsv'

    shell:
        """
        wget '{params.gex}' -O '{output.gex}'
        wget '{params.peaks}' -O '{output.peaks}'
        wget '{params.frags}' -O '{output.frags}'
        """

rule index_frags:
    input:
        frags='datasets/pitunpaired/smpl.frags.tsv.gz'
    output:
        index='datasets/pitunpaired/smpl.frags.tsv.gz.tbi'
    params:
        unzip='datasets/pitunpaired/smpl.frags.tsv'

    singularity:
        'workflow/envs/gretabench.sif'
    
    shell:
        """
        gunzip -d {input.frags}
        bgzip {params.unzip}
        tabix -p bed {input.frags}
        """


rule coembedd_pitunpaired:
    input:
        gex='datasets/pitunpaired/smpl.filtered_feature_bc_matrix.h5',
        celltypes='workflow/scripts/datasets/pitunpaired/celltypes.csv',
        peaks='datasets/pitunpaired/peaks.original.h5',
        frags='datasets/pitunpaired/smpl.frags.tsv.gz'

    
    output:
        annot=temp(local('datasets/pitunpaired/annot.csv'))
        exprMat=temp(local('datasets/pitunpaired/exprMat.rds')),
        atacSE=temp(local('datasets/pitunpaired/atac.se.rds')),
        cca=temp(local('datasets/pitunpaired/cca.rds'))

    singularity: 
        'workflow/envs/seurat.sif'

    shell:
        """
        Rscript workflow/scripts/datasets/pitunpaired/coembedd.R \
        {input.gex} \
        {input.celltypes} \
        {input.peaks} \
        {input.frags} \
        {output.annot} \
        {output.exprMat} \
        {output.atacSE} \
        {output.cca} 
        """


rule pairCells_pitunpaired:
    input:
        exprMat=local('datasets/pitunpaired/exprMat.rds'),
        atacSE=local('datasets/pitunpaired/atac.se.rds'),
        cca=local('datasets/pitunpaired/cca.rds')
    output:
        barMap=local('datasets/pitunpaired/barMap.csv')

    singularity:
        'workflow/envs/figr.sif'

    shell:
        """
        Rscript workflow/scripts/datasets/pitunpaired/pairCells.R \
        {input.exprMat} \
        {input.atacSE} \
        {input.cca} \
        {output.barMap} \
        """

rule callpeaks_pitunpaired:
    input:
        frags='datasets/pitunpaired/smpl.frags.tsv.gz',
        annot='datasets/pitunpaired/annot.csv',
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        tmp=temp(directory(local('datasets/pitunpaired/tmp_peaks'))),
        peaks=local('datasets/pitunpaired/peaks.h5ad')
    resources:
        mem_mb=32000,
    threads: 16
    shell:
        """
        python workflow/scripts/datasets/callpeaks.py \
        -f {input.frags} \
        -a {input.annot} \
        -t {output.tmp} \
        -o {output.peaks}
        """

rule annotate_pitunpaired:
    input:
        annot='workflow/scripts/datasets/pitunpaired/annot.csv',
        g='gdata/geneids',
        peaks='datasets/pitunpaired/peaks.h5ad',
        gex='datasets/pitunpaired/smpl.filtered_feature_bc_matrix.h5',
        barmap='datasets/pitunpaired/barMap.csv'
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        tmp=temp(directory(local('datasets/pitunpaired/tmp_annot'))),
        out='datasets/pitunpaired/annotated.h5mu'
    params:
        organism=config['datasets']['pitunpaired']['organism'],
    resources:
        mem_mb=32000,
    shell:
        """
        python workflow/scripts/datasets/pitunpaired/pitunpaired.py \
        -a {output.tmp} \
        -b {input.annot} \
        -c {input.g} \
        -d {params.organism} \
        -e {input.peaks} \
        -f {output.out} \
        -g {input.gex} \
        -i {input.barmap}
        """

