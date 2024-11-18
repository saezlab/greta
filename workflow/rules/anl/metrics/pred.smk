rule pred_omics:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        grn=lambda w: rules.grn_run.output.out.format(**w),
    output:
        out='anl/metrics/pred/omics/{db}/{dat}.{case}/{pre}.{p2g}.{tfb}.{mdl}.scores.csv'
    params:
        col_source=lambda w: 'cre' if w.db == 'gcre' else 'source',
        col_target=lambda w: 'cre' if w.db == 'cretf' else 'target',
        mod_source=lambda w: 'atac' if w.db == 'gcre' else 'rna',
        mod_target=lambda w: 'atac' if w.db == 'cretf' else 'rna',
    shell:
        """
        python workflow/scripts/anl/metrics/pred/omics.py \
        -a {input.grn} \
        -b {params.col_source} \
        -c {params.col_target} \
        -d {params.mod_source} \
        -e {params.mod_target} \
        -f {output}
        """


rule pred_pathway:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        grn=lambda w: rules.grn_run.output.out.format(**w),
        rsc=lambda w: expand('dbs/hg38/gst/{gst_db}.csv', gst_db=config['dbs'][w.org]['gst'][w.gst_db])
    output:
        out='analysis/metrics/pred/gsets/{gst_db}/{dat}.{case}/{pre}.{p2g}.{tfb}.{mdl}.scores.csv'
    shell:
        """
        python workflow/scripts/anl/metrics/pred/gsets.py \
        -i {input.grn} \
        -p {input.rsc} \
        -o {output}
        """
