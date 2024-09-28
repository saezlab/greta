localrules: download_pitunpair, index_frags


rule download_pitunpair:
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
        wget '{params.gex}' -O '{output.gex}'
        wget '{params.peaks}' -O '{output.peaks}'
        wget '{params.frags}' -O '{output.frags}'
        wget '{params.celltypes}' -O '{output.celltypes}'
        """


rule index_frags:
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
        exprmat=temp(local('datasets/pitunpair/exprmat.rds')),
        atacse=temp(local('datasets/pitunpair/atacse.rds')),
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
        {output.exprmat} \
        {output.atacse} \
        {output.cca}
        """


rule paircells_pitunpair:
    threads: 32
    input:
        exprmat=rules.coembedd_pitunpair.output.exprmat,
        atacse=rules.coembedd_pitunpair.output.atacse,
        cca=rules.coembedd_pitunpair.output.cca,
        celltypes=rules.download_pitunpair.output.celltypes,
    output:
        barmap='datasets/pitunpair/barmap.csv'
    singularity:
        'workflow/envs/figr.sif'
    shell:
        """
        Rscript workflow/scripts/datasets/pitunpair/pairCells.R \
        {input.exprMat} \
        {input.atacSE} \
        {input.cca} \
        {input.celltypes} \
        {output.barMap} \
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
        -o {output.peaks}
        """


rule annotate_pitunpair:
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

