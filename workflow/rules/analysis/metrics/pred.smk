rule pred_omics_gtf:
    input:
        'datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.{mdl}.grn.csv'
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        'analysis/metrics/pred/omics/gtf/{dataset}.{case}/{pre}.{p2g}.{tfb}.{mdl}.scores.csv'
    params:
        col_source='source',
        col_target='target',
        mod_source='rna',
        mod_target='rna',
    shell:
        """
        python workflow/scripts/analysis/metrics/pred/compute_omics.py \
        -a {input} \
        -b {params.col_source} \
        -c {params.col_target} \
        -d {params.mod_source} \
        -e {params.mod_target} \
        -f {output}
        """


rule pred_omics_cretf:
    input:
        'datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.{mdl}.grn.csv'
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        'analysis/metrics/pred/omics/cretf/{dataset}.{case}/{pre}.{p2g}.{tfb}.{mdl}.scores.csv'
    params:
        col_source='source',
        col_target='cre',
        mod_source='rna',
        mod_target='atac',
    shell:
        """
        python workflow/scripts/analysis/metrics/pred/compute_omics.py \
        -a {input} \
        -b {params.col_source} \
        -c {params.col_target} \
        -d {params.mod_source} \
        -e {params.mod_target} \
        -f {output}
        """


rule pred_omics_gcre:
    input:
        'datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.{mdl}.grn.csv'
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        'analysis/metrics/pred/omics/gcre/{dataset}.{case}/{pre}.{p2g}.{tfb}.{mdl}.scores.csv'
    params:
        col_source='cre',
        col_target='target',
        mod_source='atac',
        mod_target='rna',
    shell:
        """
        python workflow/scripts/analysis/metrics/pred/compute_omics.py \
        -a {input} \
        -b {params.col_source} \
        -c {params.col_target} \
        -d {params.mod_source} \
        -e {params.mod_target} \
        -f {output}
        """


rule pred_pathway_kegg:
    input:
        'datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.{mdl}.grn.csv'
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        'analysis/metrics/pred/pathway/kegg/{dataset}.{case}/{pre}.{p2g}.{tfb}.{mdl}.scores.csv'
    params:
        kegg_path='/mnt/sds-hd/sd22b002/projects/GRETA/greta_resources/database/hg38/spred/pways/allkegg.csv'
    shell:
        """
        python workflow/scripts/analysis/metrics/pred/compute_ptwpred.py \
        -i {input} \
        -p {params.kegg_path} \
        -o {output}
        """


rule pred_pathway_hall:
    input:
        'datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.{mdl}.grn.csv'
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        'analysis/metrics/pred/pathway/hall/{dataset}.{case}/{pre}.{p2g}.{tfb}.{mdl}.scores.csv'
    params:
        hall_path='/mnt/sds-hd/sd22b002/projects/GRETA/greta_resources/database/hg38/spred/pways/hall.csv'
    shell:
        """
        python workflow/scripts/analysis/metrics/pred/compute_ptwpred.py \
        -i {input} \
        -p {params.hall_path} \
        -o {output}
        """


rule pred_pathway_reac:
    input:
        'datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.{mdl}.grn.csv'
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        'analysis/metrics/pred/pathway/reac/{dataset}.{case}/{pre}.{p2g}.{tfb}.{mdl}.scores.csv'
    params:
        reac_path='/mnt/sds-hd/sd22b002/projects/GRETA/greta_resources/database/hg38/spred/pways/reac.csv'
    shell:
        """
        python workflow/scripts/analysis/metrics/pred/compute_ptwpred.py \
        -i {input} \
        -p {params.reac_path} \
        -o {output}
        """
