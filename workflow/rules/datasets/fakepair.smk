rule index_frags_fakepair:
    threads: 1
    input:
        frags=lambda w: map_rules('download_', rule_name='{dname}pair', w.dname, out='frags'),
    output:
        frags=temp(local('datasets/fake{dname}pair/smpl.frags.tsv.gz')),
        index=temp(local('datasets/fake{dname}pair/smpl.frags.tsv.gz.tbi')),
    params:
        unzip='datasets/fake{dname}pair/smpl.frags.tsv'
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
        gex=lambda w: map_rules('download_', rule_name='{dname}pair', w.dname, out='gex'),
        peaks=lambda w: map_rules('callpeaks_', rule_name='{dname}pair', w.dname, out='peaks'),
        frags=rules.index_frags_fakepair.output.frags,
        index=rules.index_frags_fakepair.output.index,
    output:
        cca=temp(local('datasets/fake{dname}pair/cca.rds'))
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
        annot=lambda w: map_rules('download_', rule_name='{dname}pair', w.dname, out='annot'),
    output:
        barmap=temp(local('datasets/fake{dname}pair/barmap.csv'))
    singularity:
        'workflow/envs/figr.sif'
    shell:
        """
        Rscript workflow/scripts/datasets/fakepair/paircells.R \
        {input.cca} \
        {input.annot} \
        {output.barmap}
        """
