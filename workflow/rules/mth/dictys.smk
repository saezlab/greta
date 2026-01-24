rule pre_dictys:
    threads: 1
    conda: '../../envs/dictys.yaml'
    input:
        mdata=rules.extract_case.output.mdata,
    output:
        out='dts/{org}/{dat}/cases/{case}/runs/dictys.pre.h5mu',
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    shell:
        """
        path_expr=$(dirname {output.out})/dictys_pre_expr.tsv.gz
        python workflow/scripts/mth/dictys/pre.py \
        -m {input.mdata} \
        -t $path_expr \
        -o {output.out} && \
        rm -f $path_expr
        """


rule p2g_dictys:
    threads: 1
    conda: '../../envs/dictys.yaml'
    input:
        pre=lambda wildcards: map_rules('pre', wildcards.pre),
        ann=lambda w: rules.gen_ann_dictys_mm10.output if config['dts'][w.dat]['organism'] == 'mm10' else rules.gen_ann_dictys.output,
    output:
        out='dts/{org}/{dat}/cases/{case}/runs/{pre}.dictys.p2g.csv',
    params:
        ext=config['methods']['dictys']['ext'] // 2,
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    shell:
        """
        path_tmp=$(dirname {output.out})/tmp_dictys_pre-{wildcards.pre}
        mkdir -p $path_tmp
        set +e
        timeout $(({resources.runtime}-20))m bash -c \
        "python workflow/scripts/mth/dictys/p2g.py \
        -d {input.pre} \
        -t $path_tmp \
        -p {output.out} \
        -e {params.ext} \
        -g {input.ann}"
        rm -rf $path_tmp
        if [ $? -eq 124 ]; then
            awk 'BEGIN {{ print "cre,gene,score" }}' > {output.out}
            rm -rf $path_tmp
        fi
        """


rule tfb_dictys:
    threads: 16
    conda: '../../envs/dictys.yaml'
    container: None
    input:
        pre=lambda wildcards: map_rules('pre', wildcards.pre),
        p2g=lambda wildcards: map_rules('p2g', wildcards.p2g),
        frags=list_frags_files,
        motif=lambda w: rules.gen_motif_dictys_mm10.output if config['dts'][w.dat]['organism'] == 'mm10' else rules.gen_motif_dictys.output,
        genome=lambda w: rules.gen_genome_dictys_mm10.output if config['dts'][w.dat]['organism'] == 'mm10' else rules.gen_genome_dictys.output,
    output:
        out='dts/{org}/{dat}/cases/{case}/runs/{pre}.{p2g}.dictys.tfb.csv'
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    params:
        use_p2g=False,
    shell:
        """
        set +e
        path_tmp=$(dirname {output.out})/tmp_dictys_pre-{wildcards.pre}_p2g-{wildcards.p2g}
        mkdir -p $path_tmp
        timeout $(({resources.runtime}-20))m \
        bash workflow/scripts/mth/dictys/tfb.sh \
        --input_pre {input.pre} \
        --input_p2g {input.p2g} \
        --input_frags {input.frags} \
        --input_motif {input.motif} \
        --input_genome {input.genome} \
        --output_d $path_tmp \
        --output_out {output.out} \
        --threads {threads} \
        --use_p2g {params.use_p2g}
        rm -rf $path_tmp
        if [ $? -eq 124 ]; then
            awk 'BEGIN {{ print "cre,tf,score" }}' > {output.out}
            rm -rf $path_tmp
        fi
        """


rule mdl_dictys:
    threads: 4
    conda: '../../envs/dictys.yaml'
    container: None
    input:
        pre=lambda wildcards: map_rules('pre', wildcards.pre),
        p2g=lambda wildcards: map_rules('p2g', wildcards.p2g),
        tfb=lambda wildcards: map_rules('tfb', wildcards.tfb),
        ann=lambda w: rules.gen_ann_dictys_mm10.output if config['dts'][w.dat]['organism'] == 'mm10' else rules.gen_ann_dictys.output,
    output:
        out='dts/{org}/{dat}/cases/{case}/runs/{pre}.{p2g}.{tfb}.dictys.mdl.csv'
    params:
        ext=config['methods']['dictys']['ext'] // 2,
        n_p2g_links=config['methods']['dictys']['n_p2g_links'],
        device=config['methods']['dictys']['device'],
        thr_score=config['methods']['dictys']['thr_score'],
        use_p2g=True,
    resources:
        partition=lambda w: "cpu-single" if w.pre=='granie' else "gpu-single",
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
        slurm=lambda w: "gres=gpu:0" if w.pre=='granie' else "gres=gpu:1",
    shell:
        """
        set +e
        path_tmp=$(dirname {output.out})/tmp_dictys_pre-{wildcards.pre}_p2g-{wildcards.p2g}_tfb-{wildcards.tfb}
        path_mdl=$path_tmp/mdl.csv
        mkdir -p $path_tmp
        timeout $(({resources.runtime}-20))m \
        bash workflow/scripts/mth/dictys/mdl_d.sh \
        --output_d $path_tmp \
        --pre_path {input.pre} \
        --p2g_path {input.p2g} \
        --tfb_path {input.tfb} \
        --mdl_path $path_mdl \
        --annot {input.ann} \
        --distance {params.ext} \
        --n_p2g_links {params.n_p2g_links} \
        --threads {threads} \
        --device {params.device} \
        --thr_score {params.thr_score} \
        --use_p2g {params.use_p2g} \
        --out_path {output.out} && \
        rm -rf $path_tmp
        if [ $? -eq 124 ]; then
            awk 'BEGIN {{ print "source,target,score,pval" }}' > {output.out}
            rm -rf $path_tmp
        fi
        """


rule mdl_o_dictys:
    threads: 4
    conda: '../../envs/dictys.yaml'
    container: None
    input:
        mdata=rules.extract_case.output.mdata,
        ann=lambda w: rules.gen_ann_dictys_mm10.output if config['dts'][w.dat]['organism'] == 'mm10' else rules.gen_ann_dictys.output,
        frags=list_frags_files,
        motif=lambda w: rules.gen_motif_dictys_mm10.output if config['dts'][w.dat]['organism'] == 'mm10' else rules.gen_motif_dictys.output,
        genome=lambda w: rules.gen_genome_dictys_mm10.output if config['dts'][w.dat]['organism'] == 'mm10' else rules.gen_genome_dictys.output,
    output:
        out='dts/{org}/{dat}/cases/{case}/runs/o_dictys.o_dictys.o_dictys.o_dictys.mdl.csv',
    params:
        ext=config['methods']['dictys']['ext'] // 2,
        n_p2g_links=config['methods']['dictys']['n_p2g_links'],
        device=config['methods']['dictys']['device'],
        thr_score=config['methods']['dictys']['thr_score'],
        use_p2g=False,
    resources:
        partition='gpu-single',
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'] * 4,
        slurm="gres=gpu:1",
    shell:
        """
        path_tmp=$(dirname {output.out})/tmp_o_dictys
        path_expr=$path_tmp/expr.tsv.gz
        path_pre=$path_tmp/pre.h5mu
        path_p2g=$path_tmp/p2g.csv
        path_tfb=$path_tmp/tfb.csv
        path_mdl=$path_tmp/mdl.csv
        set +e
        timeout $(({resources.runtime}-20))m bash -c \
        "mkdir -p $path_tmp && \
        python workflow/scripts/mth/dictys/pre.py \
        -m {input.mdata} \
        -t $path_expr \
        -o $path_pre && \
        python workflow/scripts/mth/dictys/p2g.py \
        -d $path_pre \
        -t $path_tmp \
        -p $path_p2g \
        -e {params.ext} \#snakemake --profile config/slurm/ --retries 0 plt/stab/fig.pdf

        -g {input.ann} && \
        bash workflow/scripts/mth/dictys/tfb.sh \
        --input_pre $path_pre \
        --input_p2g $path_p2g \
        --input_frags {input.frags} \
        --input_motif {input.motif} \
        --input_genome {input.genome} \
        --output_d $path_tmp \
        --output_out $path_tfb \
        --use_p2g {params.use_p2g} \
        --threads {threads} && \
        bash workflow/scripts/mth/dictys/mdl.sh \
        --output_d $path_tmp \
        --pre_path $path_pre \
        --p2g_path $path_p2g \
        --tfb_path $path_tfb \
        --mdl_path $path_mdl \
        --annot {input.ann} \
        --distance {params.ext} \
        --n_p2g_links {params.n_p2g_links} \
        --threads {threads} \
        --device {params.device} \
        --thr_score {params.thr_score} \
        --use_p2g {params.use_p2g} \
        --out_path {output.out} && \
        rm -rf $path_tmp"
        if [ $? -eq 124 ]; then
            awk 'BEGIN {{ print "source,cre,target,score,pval" }}' > {output.out}
            rm -rf $path_tmp
        fi
        """
