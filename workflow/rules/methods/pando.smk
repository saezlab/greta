localrules: download_granges


rule download_granges:
    output:
        h='gdata/granges/hg38_ensdb_v86.csv',
        m='gdata/granges/mm10_ensdb_v79.csv',
        d=directory('gdata/granges'),
    singularity:
        'workflow/envs/pando.sif'
    shell:
        """
        Rscript workflow/scripts/methods/pando/get_granges.R {output.h} {output.m}
        """


rule pre_pando:
    threads: 32
    input:
        mdata=rules.extract_case.output.mdata,
        h=rules.download_granges.output.h,
        m=rules.download_granges.output.m,
    singularity:
        'workflow/envs/pando.sif'
    output:
        p=temp(local('datasets/{dataset}/cases/{case}/runs/pando.peaks.csv')),
        m=temp(local('datasets/{dataset}/cases/{case}/runs/pando.matches.csv')),
        out='datasets/{dataset}/cases/{case}/runs/pando.pre.h5mu'
    params:
        organism=lambda w: config['datasets'][w.dataset]['organism'],
        exclude_exons=config['methods']['pando']['exclude_exons'],
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    shell:
        """
        Rscript workflow/scripts/methods/pando/pre.R \
        {input.mdata} \
        {params.organism} \
        {input.h} \
        {input.m} \
        {params.exclude_exons} \
        {output.p} \
        {output.m}
        python workflow/scripts/methods/pando/pre.py \
        -i {input.mdata} \
        -p {output.p} \
        -m {output.m} \
        -o {output.out}
        """


rule p2g_pando:
    threads: 32
    input:
        pre=lambda wildcards: map_rules('pre', wildcards.pre),
        h=rules.download_granges.output.h,
        m=rules.download_granges.output.m,
    singularity:
        'workflow/envs/pando.sif'
    output:
        out='datasets/{dataset}/cases/{case}/runs/{pre}.pando.p2g.csv'
    params:
        organism=lambda w: config['datasets'][w.dataset]['organism'],
        ext=config['methods']['pando']['ext'],
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    shell:
        """
        set +e
        timeout $(({resources.runtime}-20))m \
        Rscript workflow/scripts/methods/pando/p2g.R \
        {input.pre} \
        {params.organism} \
        {input.h} \
        {input.m} \
        {params.ext} \
        {output.out}
        if [ $? -eq 124 ]; then
            awk 'BEGIN {{ print "cre,gene,score,pval" }}' > {output.out}
        fi
        """


rule tfb_pando:
    threads: 32
    input:
        pre=lambda wildcards: map_rules('pre', wildcards.pre),
        p2g=lambda wildcards: map_rules('p2g', wildcards.p2g),
    singularity:
        'workflow/envs/pando.sif'
    output:
        out='datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.pando.tfb.csv'
    params:
        organism=lambda w: config['datasets'][w.dataset]['organism']
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    shell:
        """
        set +e
        timeout $(({resources.runtime}-20))m \
        Rscript workflow/scripts/methods/pando/tfb.R \
        {input.pre} \
        {params.organism} \
        {input.p2g} \
        {output.out}
        if [ $? -eq 124 ]; then
            awk 'BEGIN {{ print "cre,tf,score" }}' > {output.out}
        fi
        """


rule mdl_pando:
    threads: 1
    input:
        pre=lambda wildcards: map_rules('pre', wildcards.pre),
        p2g=lambda wildcards: map_rules('p2g', wildcards.p2g),
        tfb=lambda wildcards: map_rules('tfb', wildcards.tfb),
    singularity:
        'workflow/envs/pando.sif'
    output:
        out='datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.pando.mdl.csv'
    params:
        thr_corr=config['methods']['pando']['thr_corr'],
        p_thresh=config['methods']['pando']['p_thresh'],
        rsq_thresh=config['methods']['pando']['rsq_thresh'],
        nvar_thresh=config['methods']['pando']['nvar_thresh'],
        min_genes_per_module=config['methods']['pando']['min_genes_per_module'],
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    shell:
        """
        set +e
        timeout $(({resources.runtime}-20))m \
        Rscript workflow/scripts/methods/pando/mdl.R \
        {input.pre} \
        {input.p2g} \
        {input.tfb} \
        {params.thr_corr} \
        {params.p_thresh} \
        {params.rsq_thresh} \
        {params.nvar_thresh} \
        {params.min_genes_per_module} \
        {threads} \
        {output.out}
        if [ $? -eq 124 ]; then
            awk 'BEGIN {{ print "source,target,score,pval" }}' > {output.out}
        fi
        """

rule mdl_o_pando:
    threads: 1
    input:
        mdata=rules.extract_case.output.mdata,
        h=rules.download_granges.output.h,
        m=rules.download_granges.output.m,
    singularity:
        'workflow/envs/pando.sif'
    output:
        out='datasets/{dataset}/cases/{case}/runs/o_pando.o_pando.o_pando.o_pando.grn.csv'
    params:
        organism=lambda w: config['datasets'][w.dataset]['organism'],
        exclude_exons=config['methods']['pando']['exclude_exons'],
        ext=config['methods']['pando']['ext'],
        thr_corr=config['methods']['pando']['thr_corr'],
        p_thresh=config['methods']['pando']['p_thresh'],
        rsq_thresh=config['methods']['pando']['rsq_thresh'],
        nvar_thresh=config['methods']['pando']['nvar_thresh'],
        min_genes_per_module=config['methods']['pando']['min_genes_per_module'],
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    shell:
        """
        set +e
        timeout $(({resources.runtime}-20))m \
        Rscript workflow/scripts/methods/pando/src.R \
        {input.mdata} \
        {params.exclude_exons} \
        {params.ext} \
        {params.thr_corr} \
        {params.p_thresh} \
        {params.rsq_thresh} \
        {params.nvar_thresh} \
        {params.min_genes_per_module} \
        {params.organism} \
        {input.h} \
        {input.m} \
        {threads} \
        {output.out}
        if [ $? -eq 124 ]; then
            awk 'BEGIN {{ print "source,target,score,pval" }}' > {output.out}
        fi
        """
