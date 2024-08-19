rule download_pitunpair:
    output:
        gex='datasets/pituunpaired/smpl.filtered_feature_bc_matrix.h5',
        peaks='datasets/pituunpaired/peaks.original.h5',
        frags='datasets/pituunpaired/smpl.frags.tsv.gz',
        fragIndex='datasets/pituunpaired/smpl.frags.tsv.gz.tbi'
    params:
        gex=config['datasets']['pituunpaired']['url']['rna_mtx'],
        peaks=config['datasets']['pituunpaired']['url']['peaks'],
        frags=config['datasets']['pituunpaired']['url']['atac_frags'],
        unzip='datasets/pituunpaired/smpl.frags.tsv'
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
        gex='datasets/pituunpaired/smpl.filtered_feature_bc_matrix.h5',
        peaks='datasets/pituunpaired/peaks.original.h5',
        frags='datasets/pituunpaired/smpl.frags.tsv.gz'
    output:
        annot=temp(local('datasets/pituunpaired/annot.csv')),
        exprMat=temp(local('datasets/pituunpaired/exprMat.rds')),
        atacSE=temp(local('datasets/pituunpaired/atac.se.rds')),
        cca=temp(local('datasets/pituunpaired/cca.rds'))
    singularity: 
        'workflow/envs/seurat.sif'
    shell:
        """
        Rscript workflow/scripts/datasets/pituunpaired/coembedd.R \
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
        exprMat=temp(local('datasets/pituunpaired/exprMat.rds')),
        atacSE=temp(local('datasets/pituunpaired/atac.se.rds')),
        cca=temp(local('datasets/pituunpaired/cca.rds'))
    output:
        barMap=local('datasets/pituunpaired/barMap.csv')
    singularity:
        'workflow/envs/figr.sif'
    shell:
        """
        Rscript workflow/scripts/datasets/pituunpaired/pairCells.R \
        {input.exprMat} \
        {input.atacSE} \
        {input.cca} \
        {output.barMap} \
        """


rule callpeaks_pitunpair:
    input:
        frags='datasets/pituunpaired/smpl.frags.tsv.gz',
        annot='datasets/pituunpaired/annot.csv',
    output:
        tmp=temp(directory(local('datasets/pituunpaired/tmp_peaks'))),
        peaks=local('datasets/pituunpaired/peaks.h5ad')
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


rule annotate_pitunpair:
    input:
        annot='datasets/pituunpaired/annot.csv',
        g='gdata/geneids',
        peaks='datasets/pituunpaired/peaks.h5ad',
        gex='datasets/pituunpaired/smpl.filtered_feature_bc_matrix.h5',
        barmap='datasets/pituunpaired/barMap.csv'
    output:
        tmp=temp(directory(local('datasets/pituunpaired/tmp_annot'))),
        out='datasets/pituunpaired/annotated.h5mu'
    params:
        organism=config['datasets']['pitunpair']['organism'],
    resources:
        mem_mb=32000,
    shell:
        """
        python workflow/scripts/datasets/pituunpaired/pituunpaired.py \
        -a {output.tmp} \
        -b {input.annot} \
        -c {input.g} \
        -d {params.organism} \
        -e {input.peaks} \
        -f {output.out} \
        -g {input.gex} \
        -i {input.barmap}
        """
