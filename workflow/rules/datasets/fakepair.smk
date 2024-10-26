rule index_frags_fakepair:
    threads: 1
    input:
        frags=rules.download_pitupair.output.frags
    output:
        frags=temp(local('datasets/fakepair/smpl.frags.tsv.gz')),
        index=temp(local('datasets/fakepair/smpl.frags.tsv.gz.tbi')),
    params:
        unzip='datasets/fakepair/smpl.frags.tsv'
    singularity:
        'workflow/envs/figr.sif'
    shell:
        """
        cp {input.frags} {output.frags}
        gunzip -d {output.frags}
        bgzip {params.unzip}
        tabix -p bed {output.frags}
        """


rule coem_fakepair:
    threads: 32
    input:
        gex=rules.download_pitupair.output.gex,
        peaks=rules.callpeaks_pitupair.output.peaks,
        frags=rules.index_frags_fakepair.output.frags,
        index=rules.index_frags_fakepair.output.index,
    output:
        cca=temp(local('datasets/fakepair/cca.rds'))
    singularity:
        'workflow/envs/figr.sif'
    resources:
        mem_mb=128000,
    shell:
        """
        Rscript workflow/scripts/datasets/fakepair/coembedd.R \
        {input.gex} \
        {input.peaks} \
        {input.frags} \
        {output.cca}
        """


rule pair_fakepair:
    threads: 1
    input:
        cca=rules.coem_fakepair.output.cca,
        annot=rules.download_pitupair.output.annot,
    output:
        barmap=temp(local('datasets/fakepair/barmap.csv'))
    singularity:
        'workflow/envs/figr.sif'
    shell:
        """
        Rscript workflow/scripts/datasets/fakepair/paircells.R \
        {input.cca} \
        {input.annot} \
        {output.barmap}
        """

rule annotate_fakepair:
    threads: 1
    input:
        annot=rules.pair_fakepair.output.barmap,
        peaks=rules.callpeaks_pitupair.output.peaks,
        gex=rules.download_pitupair.output.gex,
        g=rules.download_geneids.output.dr,
    output:
        out='datasets/fakepair/annotated.h5mu'
    singularity:
        'workflow/envs/gretabench.sif'
    params:
        organism=config['datasets']['pitupair']['organism'],
    resources:
        mem_mb=32000,
    shell:
        """
        python workflow/scripts/datasets/fakepair/fakepair.py \
        -b {input.annot} \
        -c {input.g} \
        -d {params.organism} \
        -e {input.peaks} \
        -f {output} \
        -g {input.gex}
        """
