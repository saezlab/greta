localrules: index_frags_fakepair


rule index_frags_fakepair:
    threads: 1
    input:
        frags=lambda w: map_rules('download', w_name='{dname}pair'.format(dname=w.dname), out='frags'),
        tbis=lambda w: map_rules('download', w_name='{dname}pair'.format(dname=w.dname), out='tbis'),
    output:
        frags=temp(local('dts/fake{dname}pair/smpl.frags.tsv.gz')),
        tbis=temp(local('dts/fake{dname}pair/smpl.frags.tsv.gz.tbi')),
    shell:
        """
        cp {input.frags} {output.frags}
        cp {input.tbis} {output.tbis}
        """


rule coem_fakepair:
    threads: 32
    singularity: 'workflow/envs/figr.sif'
    input:
        gex=lambda w: map_rules(rule_prefix='download', w_name='{dname}pair'.format(dname=w.dname), out='gex'),
        peaks=lambda w: map_rules('callpeaks', w_name='{dname}pair'.format(dname=w.dname), out='peaks'),
        frags=rules.index_frags_fakepair.output.frags,
        tbis=rules.index_frags_fakepair.output.tbis,
    output:
        cca=temp(local('dts/fake{dname}pair/cca.rds'))
    resources: mem_mb=128000
    shell:
        """
        Rscript workflow/scripts/dts/fakepair/coembedd.R \
        {input.gex} \
        {input.peaks} \
        {input.frags} \
        {output.cca}
        """


rule pair_fakepair:
    threads: 1
    singularity: 'workflow/envs/figr.sif'
    input:
        cca=rules.coem_fakepair.output.cca,
        annot=lambda w: map_rules(rule_prefix='download', w_name='{dname}pair'.format(dname=w.dname), out='annot'),
    output: barmap=temp(local('dts/fake{dname}pair/barmap.csv'))
    shell:
        """
        Rscript workflow/scripts/dts/fakepair/paircells.R \
        {input.cca} \
        {input.annot} \
        {output.barmap}
        """

localrules: annotate_fakepitupair
rule annotate_fakepitupair:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        mdata=rules.annotate_pitupair.output.out,
        barmap='dts/fakepitupair/barmap.csv',
    output:
        out='dts/fakepitupair/annotated.h5mu'
    shell:
        """
        python workflow/scripts/dts/fakepair/fakepair.py \
        -m {input.mdata} \
        -b {input.barmap} \
        -o {output.out}
        """
