rule pred_omics:
    input:
        grn=lambda wildcards: rules.grn_run.output.out.format(**wildcards),
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        out='analysis/metrics/pred/omics/{resource}/{dataset}.{case}/{pre}.{p2g}.{tfb}.{mdl}.scores.csv'
    params:
        col_source=lambda w: 'cre' if w.resource == 'gcre' else 'source',
        col_target=lambda w: 'cre' if w.resource == 'cretf' else 'target',
        mod_source=lambda w: 'atac' if w.resource == 'gcre' else 'rna',
        mod_target=lambda w: 'atac' if w.resource == 'cretf' else 'rna',
    shell:
        """
        python workflow/scripts/analysis/metrics/pred/compute_omics.py \
        -a {input.grn} \
        -b {params.col_source} \
        -c {params.col_target} \
        -d {params.mod_source} \
        -e {params.mod_target} \
        -f {output}
        """


rule pred_pathway:
    input:
        grn=lambda wildcards: rules.grn_run.output.out.format(**wildcards),
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        out='analysis/metrics/pred/pathway/{resource}/{dataset}.{case}/{pre}.{p2g}.{tfb}.{mdl}.scores.csv'
    params:
        pw_path='/mnt/sds-hd/sd22b002/projects/GRETA/greta_resources/database/hg38/spred/pways/{resource}.csv'
    shell:
        """
        python workflow/scripts/analysis/metrics/pred/compute_ptwpred.py \
        -i {input.grn} \
        -p {params.pw_path} \
        -o {output}
        """


rule pred_pair_pitu:
    input:
        p='analysis/topo/{dname}pair.{case}.sims_mult.csv',
        n='analysis/topo/{dname}npair.{case}.sims_mult.csv',
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        make_combs(
            path='analysis/metrics/pred/pair/{dname}/{dname}.{case}/',
            mthds=mthds,
            name='scores',
        )
    shell:
        """
        python workflow/scripts/analysis/metrics/pred/compute_pair.py -i {input.p}
        """

