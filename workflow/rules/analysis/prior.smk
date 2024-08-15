def make_combs(path):
    mthds = ['celloracle', 'figr', 'granie', 'pando']
    from itertools import product
    s = '{0}.{1}.{2}.{3}.scores.csv'
    combinations = product(mthds, repeat=4)
    strings = []
    for combo in combinations:
        strings.append(path + '{dataset}.{case}.' + s.format(*combo))
    s = path + '{dataset}.{case}.random.random.random.random.scores.csv'
    strings.append(s)
    return strings

rule aggr_type_task_resource:
    input:
        make_combs(path='analysis/metrics/{type}/{task}/{resource}/')
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        'analysis/metrics/{type}/{task}/{resource}/{dataset}.{case}.scores.csv'
    shell:
        """
        python workflow/scripts/analysis/metrics/aggregate.py \
        -i {input} \
        -o {output}
        """

rule prior_tfm:
    input:
        grn='datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.{mdl}.grn.csv',
        resource='/mnt/sds-hd/sd22b002/projects/GRETA/greta_resources/database/hg38/tfm/{resource}/{resource}.csv',
        cats='config/prior_cats.json',
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        temp(local('analysis/metrics/prior/tfm/{resource}/{dataset}.{case}.{pre}.{p2g}.{tfb}.{mdl}.scores.csv'))
    shell:
        """
        python workflow/scripts/analysis/metrics/prior/compute_tfm.py \
        -a {input.grn} \
        -b {input.resource} \
        -c {input.cats} \
        -f {output}
        """

rule prior_tfbind:
    input:
        grn='datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.{mdl}.grn.csv',
        resource='/mnt/sds-hd/sd22b002/projects/GRETA/greta_resources/database/hg38/tfbind/{resource}/{resource}.bed',
        cats='config/prior_cats.json',
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        temp(local('analysis/metrics/prior/tfbind/{resource}/{dataset}.{case}.{pre}.{p2g}.{tfb}.{mdl}.scores.csv'))
    params:
        grp='source',
    shell:
        """
        python workflow/scripts/analysis/metrics/prior/compute_gnm.py \
        -a {input.grn} \
        -b {input.resource} \
        -c {input.cats} \
        -d {params.grp} \
        -f {output}
        """

rule prior_cre:
    input:
        grn='datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.{mdl}.grn.csv',
        resource='/mnt/sds-hd/sd22b002/projects/GRETA/greta_resources/database/hg38/cre/{resource}/{resource}.bed',
        cats='config/prior_cats.json',
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        temp(local('analysis/metrics/prior/cre/{resource}/{dataset}.{case}.{pre}.{p2g}.{tfb}.{mdl}.scores.csv'))
    shell:
        """
        python workflow/scripts/analysis/metrics/prior/compute_gnm.py \
        -a {input.grn} \
        -b {input.resource} \
        -c {input.cats} \
        -f {output}
        """

rule prior_eqtl:
    input:
        grn='datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.{mdl}.grn.csv',
        resource='/mnt/sds-hd/sd22b002/projects/GRETA/greta_resources/database/hg38/eqtl/{resource}/{resource}.bed',
        cats='config/prior_cats.json',
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        temp(local('analysis/metrics/prior/eqtl/{resource}/{dataset}.{case}.{pre}.{p2g}.{tfb}.{mdl}.scores.csv'))
    params:
        grp='target',
    shell:
        """
        python workflow/scripts/analysis/metrics/prior/compute_gnm.py \
        -a {input.grn} \
        -b {input.resource} \
        -c {input.cats} \
        -d {params.grp} \
        -f {output}
        """
