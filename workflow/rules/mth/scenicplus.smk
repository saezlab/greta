rule mdl_o_scenicplus:
    threads: 32
    singularity: 'workflow/envs/scenicplus.sif'
    input:
        img='workflow/envs/scenicplus.sif',
        mdata=rules.extract_case.output.mdata,
        blist=rules.cre_blacklist.output,
        rnk=rules.gen_motif_scenicplus.output.human_rankings,
        man=rules.gen_motif_scenicplus.output.human_annot,
        scr=rules.gen_motif_scenicplus.output.human_scores,
        ann=rules.gen_genome_scenicplus.output.ann,
        csz=rules.gen_genome_scenicplus.output.csz,
    output:
        dir=directory('dts/{dat}/cases/{case}/runs/scenicplus/'),
        out='dts/{dat}/cases/{case}/runs/o_scenicplus.o_scenicplus.o_scenicplus.o_scenicplus.mdl.csv'
    params:
        ntopics=config['methods']['scenicplus']['ntopics'],
        ext=config['methods']['scenicplus']['ext'] // 2,
    resources:
        mem_mb=lambda wildcards, attempt: restart_mem(wildcards, attempt) * 4,
        runtime=config['max_mins_per_step'],
    shell:
        """
        mkdir -p {output.dir}
        set +e
        timeout $(({resources.runtime}-20))m \
        bash workflow/scripts/mth/scenicplus/o_mdl.sh \
        --new_dir {output.dir} \
        --path_mdata {input.mdata} \
        --path_blist {input.blist} \
        --ntopics {params.ntopics} \
        --path_ann {input.ann} \
        --path_csz {input.csz} \
        --ext {params.ext} \
        --path_rnk {input.rnk} \
        --path_man {input.man} \
        --path_scr {input.scr} \
        --threads {threads} \
        --path_out {output.out}
        if [ $? -eq 124 ]; then
            awk 'BEGIN {{ print "source,target,score,pval" }}' > {output.out}
        fi
        """


rule pre_scenicplus:
    threads: 1
    singularity: 'workflow/envs/scenicplus.sif'
    input:
        mdata=rules.extract_case.output.mdata,
        dir=rules.mdl_o_scenicplus.output.dir,
    output:
        out='dts/{dat}/cases/{case}/runs/scenicplus.pre.h5mu'
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    shell:
        """
        python workflow/scripts/mth/scenicplus/pre.py {input.mdata} {input.dir}/mdata.h5mu {output.out}
        """


rule p2g_scenicplus:
    threads: 32
    singularity: 'workflow/envs/scenicplus.sif'
    input:
        dir=rules.mdl_o_scenicplus.output.dir,
        pre=lambda wildcards: map_rules('pre', wildcards.pre),
        ann=rules.gen_genome_scenicplus.output.ann,
        csz=rules.gen_genome_scenicplus.output.csz,
    output:
        out='dts/{dat}/cases/{case}/runs/{pre}.scenicplus.p2g.csv'
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    params:
        ext=config['methods']['scenicplus']['ext'] // 2,
    shell:
        """
        new_dir=$(dirname {output.out})/{wildcards.pre}_scenicplus_p2g/
        mkdir -p $new_dir
        set +e
        timeout $(({resources.runtime}-20))m \
        bash workflow/scripts/mth/scenicplus/p2g.sh \
        --new_dir $new_dir \
        --path_pre {input.pre} \
        --path_ann {input.ann} \
        --path_csz {input.csz} \
        --ext {params.ext} \
        --threads {threads} \
        --path_out {output.out}
        if [ $? -eq 124 ]; then
            awk 'BEGIN {{ print "cre,gene,score" }}' > {output.out}
        fi
        rm -rf $new_dir
        """


rule tfb_scenicplus:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        dir=rules.mdl_o_scenicplus.output.dir,
        pre=lambda wildcards: map_rules('pre', wildcards.pre),
        p2g=lambda wildcards: map_rules('p2g', wildcards.p2g),
    output:
        out='dts/{dat}/cases/{case}/runs/{pre}.{p2g}.scenicplus.tfb.csv'
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    shell:
        """
        set +e
        timeout $(({resources.runtime}-20))m bash -c \
        'python workflow/scripts/mth/scenicplus/tfb.py \
        {input.pre} {input.p2g} {input.dir}/direct.h5ad {output.out}'
        if [ $? -eq 124 ]; then
            awk 'BEGIN {{ print "cre,tf,score" }}' > {output.out}
        fi
        """


rule mdl_scenicplus:
    threads: 32
    singularity: 'workflow/envs/scenicplus.sif'
    input:
        pre=lambda wildcards: map_rules('pre', wildcards.pre),
        p2g=lambda wildcards: map_rules('p2g', wildcards.p2g),
        tfb=lambda wildcards: map_rules('tfb', wildcards.tfb),
        rnk=rules.gen_motif_scenicplus.output.human_rankings,
    output:
        out='dts/{dat}/cases/{case}/runs/{pre}.{p2g}.{tfb}.scenicplus.mdl.csv'
    resources:
        mem_mb=lambda wildcards, attempt: restart_mem(wildcards, attempt) * 2,
        runtime=config['max_mins_per_step'],
    shell:
        """
        new_dir=$(dirname {output.out})/{wildcards.pre}_{wildcards.p2g}_{wildcards.tfb}_scenicplus_mdl/
        mkdir -p $new_dir
        set +e
        timeout $(({resources.runtime}-20))m \
        bash workflow/scripts/mth/scenicplus/mdl.sh \
        --new_dir $new_dir \
        --path_pre {input.pre} \
        --path_p2g {input.p2g} \
        --path_tfb {input.tfb} \
        --path_rnk {input.rnk} \
        --threads {threads} \
        --path_out {output.out} && \
        rm -rf $new_dir
        if [ $? -eq 124 ]; then
            awk 'BEGIN {{ print "source,target,score,pval" }}' > {output.out}
        fi
        
        """


