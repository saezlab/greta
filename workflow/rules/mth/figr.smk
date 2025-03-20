rule pre_figr:
    threads: 1
    singularity: 'workflow/envs/figr.sif'
    input:
        img='workflow/envs/figr.sif',
        mdata=rules.extract_case.output.mdata
    output:
        out='dts/{dat}/cases/{case}/runs/figr.pre.h5mu'
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    shell:
        """
        cp {input.mdata} {output.out}
        Rscript workflow/scripts/mth/figr/pre.R \
        {output.out} \
        {threads}
        """


rule p2g_figr:
    threads: 1
    singularity: 'workflow/envs/figr.sif'
    input:
        pre=lambda wildcards: map_rules('pre', wildcards.pre),
    output:
        out='dts/{dat}/cases/{case}/runs/{pre}.figr.p2g.csv'
    params:
        organism=lambda w: config['dts'][w.dat]['organism'],
        ext=config['methods']['figr']['ext'],
        thr_p2g_pval=config['methods']['figr']['thr_p2g_pval'],
        ncres=config['methods']['figr']['ncres'],
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    shell:
        """
        set +e
        timeout $(({resources.runtime}-20))m \
        Rscript workflow/scripts/mth/figr/p2g.R \
        {input.pre} \
        {params.organism} \
        {params.ext} \
        {params.thr_p2g_pval} \
        {params.ncres} \
        {threads} \
        {output.out}
        if [ $? -eq 124 ]; then
            awk 'BEGIN {{ print "cre,gene,score,pval" }}' > {output.out}
        fi
        """


rule tfb_figr:
    threads: 1
    singularity: 'workflow/envs/figr.sif'
    input:
        pre=lambda wildcards: map_rules('pre', wildcards.pre),
        p2g=lambda wildcards: map_rules('p2g', wildcards.p2g),
    output:
        out='dts/{dat}/cases/{case}/runs/{pre}.{p2g}.figr.tfb.csv'
    params:
        organism=lambda w: config['dts'][w.dat]['organism'],
        cellK=config['methods']['figr']['cellK'],
        dorcK=config['methods']['figr']['dorcK'],
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    shell:
        """
        set +e
        timeout $(({resources.runtime}-20))m \
        Rscript workflow/scripts/mth/figr/tfb.R \
        {input.pre} \
        {params.organism} \
        {input.p2g} \
        {params.cellK} \
        {params.dorcK} \
        {threads} \
        {output.out}
        if [ $? -eq 124 ]; then
            awk 'BEGIN {{ print "cre,tf,score" }}' > {output.out}
        fi
        """


rule mdl_figr:
    threads: 1
    singularity: 'workflow/envs/figr.sif'
    input:
        pre=lambda wildcards: map_rules('pre', wildcards.pre),
        p2g=lambda wildcards: map_rules('p2g', wildcards.p2g),
        tfb=lambda wildcards: map_rules('tfb', wildcards.tfb),
    output:
        out='dts/{dat}/cases/{case}/runs/{pre}.{p2g}.{tfb}.figr.mdl.csv'
    params:
        cellK=config['methods']['figr']['cellK'],
        thr_score=config['methods']['figr']['thr_score'],
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    shell:
        """
        set +e
        timeout $(({resources.runtime}-20))m \
        Rscript workflow/scripts/mth/figr/mdl.R \
        {input.pre} \
        {input.p2g} \
        {input.tfb} \
        {params.cellK} \
        {params.thr_score} \
        {threads} \
        {output.out}
        if [ $? -eq 124 ]; then
            awk 'BEGIN {{ print "source,target,score,pval" }}' > {output.out}
        fi
        """


rule mdl_o_figr:
    threads: 1
    singularity: 'workflow/envs/figr.sif'
    input: rules.extract_case.output.mdata,
    output:
        out='dts/{dat}/cases/{case}/runs/o_figr.o_figr.o_figr.o_figr.mdl.csv'
    params:
        organism=lambda w: config['dts'][w.dat]['organism'],
        ext=config['methods']['figr']['ext'],
        thr_p2g_pval=config['methods']['figr']['thr_p2g_pval'],
        ncres=config['methods']['figr']['ncres'],
        cellK=config['methods']['figr']['cellK'],
        dorcK=config['methods']['figr']['dorcK'],
        thr_score=config['methods']['figr']['thr_score'],
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'] * 2,
    shell:
        """
        set +e
        timeout $(({resources.runtime}-20))m \
        Rscript workflow/scripts/mth/figr/src.R \
        {input} \
        {params.cellK} \
        {params.organism} \
        {params.ext} \
        {params.thr_p2g_pval} \
        {params.ncres} \
        {params.dorcK} \
        {params.thr_score} \
        {threads} \
        {output.out}
        if [ $? -eq 124 ]; then
            awk 'BEGIN {{ print "source,target,score,pval" }}' > {output.out}
        fi
        """
