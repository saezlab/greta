localrules: install_dictys


rule pre_dictys:
    threads: 1
    conda: '{home_path}/miniforge3/envs/dictys'.format(home_path=home_path)
    input:
        img=rules.install_dictys.output,
        mdata=rules.extract_case.output.mdata
    output:
        tmp='dts/{dat}/cases/{case}/runs/dictys_pre_expr.tsv.gz',
        out='dts/{dat}/cases/{case}/runs/dictys.pre.h5mu',
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    shell:
        """
        python workflow/scripts/mth/dictys/pre.py \
        -m {input.mdata} \
        -t {output.tmp} \
        -o {output.out}
        """


rule p2g_dictys:
    threads: 1
    conda: '{home_path}/miniforge3/envs/dictys'.format(home_path=home_path)
    input:
        pre=lambda wildcards: map_rules('pre', wildcards.pre),
        ann=rules.gen_ann_dictys.output,
    output:
        d=temp(directory(local('dts/{dat}/cases/{case}/runs/{pre}_dictys_tmp/'))),
        out='dts/{dat}/cases/{case}/runs/{pre}.dictys.p2g.csv',
    params:
        ext=config['methods']['dictys']['ext'] // 2,
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    shell:
        """
        mkdir -p {output.d}
        set +e
        timeout $(({resources.runtime}-20))m bash -c \
        'python workflow/scripts/mth/dictys/p2g.py \
        -d {input.pre} \
        -t {output.d} \
        -p {output.out} \
        -e {params.ext} \
        -g {input.ann}'
        if [ $? -eq 124 ]; then
            awk 'BEGIN {{ print "cre,gene,score" }}' > {output.out}
        fi
        """


rule tfb_dictys:
    threads: 16
    conda: '{home_path}/miniforge3/envs/dictys'.format(home_path=home_path)
    container: None
    input:
        pre=lambda wildcards: map_rules('pre', wildcards.pre),
        p2g=lambda wildcards: map_rules('p2g', wildcards.p2g),
        frags=list_frags_files,
        motif=rules.gen_motif_dictys.output,
        genome=rules.gen_genome_dictys.output,
    output:
        d=temp(directory('dts/{dat}/cases/{case}/runs/{pre}.{p2g}.dictys_tmp')),
        out='dts/{dat}/cases/{case}/runs/{pre}.{p2g}.dictys.tfb.csv'
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    params:
        use_p2g=False,
    shell:
        """
        set +e
        timeout $(({resources.runtime}-20))m \
        bash workflow/scripts/mth/dictys/tfb.sh \
        --input_pre {input.pre} \
        --input_p2g {input.p2g} \
        --input_frags {input.frags} \
        --input_motif {input.motif} \
        --input_genome {input.genome} \
        --output_d {output.d} \
        --output_out {output.out} \
        --threads {threads} \
        --use_p2g {params.use_p2g} 
        if [ $? -eq 124 ]; then
            awk 'BEGIN {{ print "cre,tf,score" }}' > {output.out}
        fi
        """


rule mdl_dictys:
    threads: 4
    conda: '{home_path}/miniforge3/envs/dictys'.format(home_path=home_path)
    container: None
    input:
        pre=lambda wildcards: map_rules('pre', wildcards.pre),
        p2g=lambda wildcards: map_rules('p2g', wildcards.p2g),
        tfb=lambda wildcards: map_rules('tfb', wildcards.tfb),
        ann=rules.gen_ann_dictys.output,
    output:
        d=temp(directory('dts/{dat}/cases/{case}/runs/{pre}.{p2g}.{tfb}.dictys_tmp')),
        out='dts/{dat}/cases/{case}/runs/{pre}.{p2g}.{tfb}.dictys.mdl.csv'
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
        timeout $(({resources.runtime}-20))m \
        bash workflow/scripts/mth/dictys/mdl.sh \
        --output_d {output.d} \
        --pre_path {input.pre} \
        --p2g_path {input.p2g} \
        --tfb_path {input.tfb} \
        --annot {input.ann} \
        --distance {params.ext} \
        --n_p2g_links {params.n_p2g_links} \
        --threads {threads} \
        --device {params.device} \
        --thr_score {params.thr_score} \
        --use_p2g {params.use_p2g} \
        --out_path {output.out}
        if [ $? -eq 124 ]; then
            awk 'BEGIN {{ print "source,target,score,pval" }}' > {output.out}
        fi
        """


rule mdl_o_dictys:
    threads: 4
    conda: '{home_path}/miniforge3/envs/dictys'.format(home_path=home_path)
    container: None
    input:
        mdata=rules.extract_case.output.mdata,
        ann=rules.gen_ann_dictys.output,
        frags=list_frags_files,
        motif=rules.gen_motif_dictys.output,
        genome=rules.gen_genome_dictys.output
    output:
        tmp='dts/{dat}/cases/{case}/runs/o_dictys_pre_expr.tsv.gz',
        d=temp(directory(local('dts/{dat}/cases/{case}/runs/o_dictys_dictys_tmp/'))),
        pre=temp(local('dts/{dat}/cases/{case}/runs/o_dictys.pre.h5mu')),
        p2g=temp(local('dts/{dat}/cases/{case}/runs/o_dictys.o_dictys.p2g.csv')),
        tfb=temp(local('dts/{dat}/cases/{case}/runs/o_dictys.o_dictys.o_dictys.tfb.csv')),
        out='dts/{dat}/cases/{case}/runs/o_dictys.o_dictys.o_dictys.o_dictys.mdl.csv',
    params:
        ext=config['methods']['dictys']['ext'] // 2,
        n_p2g_links=config['methods']['dictys']['n_p2g_links'],
        device=config['methods']['dictys']['device'],
        thr_score=config['methods']['dictys']['thr_score'],
        use_p2g=False,
    resources:
        partition='gpu-single',
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'] * 2,
        slurm="gres=gpu:1",
    shell:
        """
        set +e
        timeout $(({resources.runtime}-20))m bash -c \
        'mkdir -p {output.d} && \
        python workflow/scripts/mth/dictys/pre.py \
        -m {input.mdata} \
        -t {output.tmp} \
        -o {output.pre} && \
        python workflow/scripts/mth/dictys/p2g.py \
        -d {output.pre} \
        -t {output.d} \
        -p {output.p2g} \
        -e {params.ext} \
        -g {input.ann} && \
        bash workflow/scripts/mth/dictys/tfb.sh \
        --input_pre {output.pre} \
        --input_p2g {output.p2g} \
        --input_frags {input.frags} \
        --input_motif {input.motif} \
        --input_genome {input.genome} \
        --output_d {output.d} \
        --output_out {output.tfb} \
        --use_p2g {params.use_p2g} \
        --threads {threads} && \
        bash workflow/scripts/mth/dictys/mdl.sh \
        --output_d {output.d} \
        --pre_path {output.pre} \
        --p2g_path {output.p2g} \
        --tfb_path {output.tfb} \
        --annot {input.ann} \
        --distance {params.ext} \
        --n_p2g_links {params.n_p2g_links} \
        --threads {threads} \
        --device {params.device} \
        --thr_score {params.thr_score} \
        --use_p2g {params.use_p2g} \
        --out_path {output.out}'
        if [ $? -eq 124 ]; then
            awk 'BEGIN {{ print "source,target,score,pval" }}' > {output.out}
        fi
        """
