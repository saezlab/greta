localrules: index_frags_fakepair


rule index_frags_fakepair:
    threads: 1
    input:
        frags=lambda w: map_rules('download', w_name='{dname}pair'.format(dname=w.dname), out='frags'),
        tbis=lambda w: map_rules('download', w_name='{dname}pair'.format(dname=w.dname), out='tbis'),
    output:
        frags=temp(local('dts/hg38/fake{dname}pair/smpl.frags.tsv.gz')),
        tbis=temp(local('dts/hg38/fake{dname}pair/smpl.frags.tsv.gz.tbi')),
    shell:
        """
        cp {input.frags} {output.frags}
        cp {input.tbis} {output.tbis}
        """


rule coem_fakepair:
    threads: 4
    conda: '../../envs/glue.yaml'
    input:
        gex=lambda w: map_rules(rule_prefix='download', w_name='{dname}pair'.format(dname=w.dname), out='gex'),
        acc=lambda w: map_rules('callpeaks', w_name='{dname}pair'.format(dname=w.dname), out='peaks'),
        gid=rules.gen_gid_ensmbl.output.hg38,
        ann=lambda w: map_rules(rule_prefix='download', w_name='{dname}pair'.format(dname=w.dname), out='annot'),
    output: temp(local('dts/hg38/fake{dname}pair/X_glue.csv'))
    resources:
        partition='gpu-single',
        slurm='gres=gpu:1'
    shell:
        """
        python workflow/scripts/dts/fakepair/coembed.py \
        {input.gex} \
        {input.acc} \
        {input.gid} \
        {input.ann} \
        {output}
        """


rule pair_fakepair:
    threads: 1
    singularity: 'workflow/envs/figr.sif'
    input:
        emd=rules.coem_fakepair.output,
        annot=lambda w: map_rules(rule_prefix='download', w_name='{dname}pair'.format(dname=w.dname), out='annot'),
    output: temp(local('dts/hg38/fake{dname}pair/barmap.csv'))
    shell:
        """
        Rscript workflow/scripts/dts/fakepair/paircells.R \
        {input.emd} \
        {input.annot} \
        {output}
        """

localrules: annotate_fakepitupair
rule annotate_fakepitupair:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        mdata=rules.annotate_pitupair.output.out,
        barmap='dts/hg38/fakepitupair/barmap.csv',
    output:
        out='dts/hg38/fakepitupair/annotated.h5mu'
    shell:
        """
        python workflow/scripts/dts/fakepair/fakepair.py \
        -m {input.mdata} \
        -b {input.barmap} \
        -o {output.out}
        """
