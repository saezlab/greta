rule pre_granie:
    threads: 1
    singularity: 'workflow/envs/granie.sif'
    input:
        img='workflow/envs/granie.sif',
        mdata=rules.extract_case.output.mdata
    output:
        out='dts/{dat}/cases/{case}/runs/granie.pre.h5mu'
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    shell:
        """
        python workflow/scripts/mth/granie/pre.py \
        -i {input.mdata} \
        -o {output.out} && \
        Rscript workflow/scripts/mth/granie/pre.R \
        {output.out} && \
        python workflow/scripts/mth/granie/pre_post.py \
        -i {output.out} \
        -o {output.out}.tmp && \
        mv {output.out}.tmp {output.out}
        """


rule p2g_granie:
    threads: 1
    singularity: 'workflow/envs/granie.sif'
    input:
        pre=lambda w: map_rules('pre', w.pre),
        gid=rules.gen_gid_ensmbl.output,
    output:
        t=temp(directory(local('dts/{dat}/cases/{case}/runs/{pre}.granie_tmp'))),
        out='dts/{dat}/cases/{case}/runs/{pre}.granie.p2g.csv'
    params: ext=config['methods']['granie']['ext'],
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    shell:
        """
        set +e
        echo {input.gid}
        timeout $(({resources.runtime}-20))m \
        Rscript workflow/scripts/mth/granie/p2g.R \
        {input.pre} \
        {input.gid} \
        {output.t} \
        {params.ext} \
        {output.out} \
        {threads} && \
        rm {output}.t
        if [ $? -eq 124 ]; then
            awk 'BEGIN {{ print "cre,gene,score,pval" }}' > {output.out}
            mkdir -p {output.t}
        fi
        """


rule tfb_granie:
    threads: 1
    singularity: 'workflow/envs/granie.sif'
    input:
        pre=lambda wildcards: map_rules('pre', wildcards.pre),
        p2g=lambda wildcards: map_rules('p2g', wildcards.p2g),
        gid=rules.gen_gid_ensmbl.output,
        tfb=rules.gen_motif_granie.output,
    output:
        t=temp(directory(local('dts/{dat}/cases/{case}/runs/{pre}.{p2g}.granie_tmp'))),
        out='dts/{dat}/cases/{case}/runs/{pre}.{p2g}.granie.tfb.csv'
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    shell:
        """
        set +e
        timeout $(({resources.runtime}-20))m \
        Rscript workflow/scripts/mth/granie/tfb.R \
        {input.pre} \
        {input.gid} \
        {output.t} \
        {input.tfb} \
        {input.p2g} \
        {output.out} \
        {threads} && \
        rm {output}.t
        if [ $? -eq 124 ]; then
            awk 'BEGIN {{ print "cre,tf,score" }}' > {output.out}
            mkdir -p {output.t}
        fi
        """


rule mdl_granie:
    threads: 1
    singularity: 'workflow/envs/granie.sif'
    input:
        pre=lambda wildcards: map_rules('pre', wildcards.pre),
        p2g=lambda wildcards: map_rules('p2g', wildcards.p2g),
        tfb=lambda wildcards: map_rules('tfb', wildcards.tfb),
        gid=rules.gen_gid_ensmbl.output,
    output:
        t=temp(directory(local('dts/{dat}/cases/{case}/runs/{pre}.{p2g}.{tfb}.granie_tmp'))),
        out='dts/{dat}/cases/{case}/runs/{pre}.{p2g}.{tfb}.granie.mdl.csv'
    params: thr_fdr=config['methods']['granie']['thr_fdr'],
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    shell:
        """
        set +e
        timeout $(({resources.runtime}-20))m \
        Rscript workflow/scripts/mth/granie/mdl.R \
        {input.pre} \
        {input.gid} \
        {output.t} \
        {input.p2g} \
        {input.tfb} \
        {params.thr_fdr} \
        {output.out} && \
        rm {output}.t
        if [ $? -eq 124 ]; then
            awk 'BEGIN {{ print "source,target,score,pval" }}' > {output.out}
            mkdir -p {output.t}
        fi
        """


rule mdl_o_granie:
    threads: 1
    singularity: 'workflow/envs/granie.sif'
    input:
        mdata=rules.extract_case.output.mdata,
        gid=rules.gen_gid_ensmbl.output,
        tfb=rules.gen_motif_granie.output,
    output:
        h=temp(local('dts/{dat}/cases/{case}/runs/pre.granie.src.h5mu')),
        t=temp(directory(local('dts/{dat}/cases/{case}/runs/granie_tmp.src'))),
        out='dts/{dat}/cases/{case}/runs/o_granie.o_granie.o_granie.o_granie.mdl.csv'
    params:
        ext=config['methods']['granie']['ext'],
        thr_fdr=config['methods']['granie']['thr_fdr'],
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'] * 2,
    shell:
        """
        set +e
        timeout $(({resources.runtime}-20))m bash -c \
        'python workflow/scripts/mth/granie/pre.py \
        -i {input.mdata} \
        -o {output.h} && \
        Rscript workflow/scripts/mth/granie/pre.R \
        {output.h} && \
        python workflow/scripts/mth/granie/pre_post.py \
        -i {output.h} \
        -o {output.h} && \
        Rscript workflow/scripts/mth/granie/src.R \
        {output.h} \
        {input.gid} \
        {output.t} \
        {params.ext} \
        {input.tfb} \
        {params.thr_fdr} \
        {output.out} \
        {threads} && \
        rm {output.t}'
        if [ $? -eq 124 ]; then
            awk 'BEGIN {{ print "source,target,score,pval" }}' > {output.out}
            mkdir -p {output.t}
        fi
        """
