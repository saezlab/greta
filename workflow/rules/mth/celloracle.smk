rule pre_celloracle:
    threads: 32
    singularity: 'workflow/envs/celloracle.sif'
    input:
        img='workflow/envs/celloracle.sif',
        mdata=rules.extract_case.output.mdata
    output:
        out='dts/{dat}/cases/{case}/runs/celloracle.pre.h5mu'
    params:
        k=config['methods']['celloracle']['k']
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    shell:
        """
        python workflow/scripts/mth/celloracle/pre.py \
        -i {input.mdata} \
        -k {params.k} \
        -o {output.out}
        """


rule p2g_celloracle:
    threads: 32
    singularity: 'workflow/envs/celloracle.sif'
    input:
        pre=lambda w: map_rules('pre', w.pre),
        csz=rules.gen_genome_celloracle.output
    output:
        pp=temp(local('dts/{dat}/cases/{case}/runs/{pre}.celloracle.peaks.csv')),
        pc=temp(local('dts/{dat}/cases/{case}/runs/{pre}.celloracle.conns.csv')),
        out='dts/{dat}/cases/{case}/runs/{pre}.celloracle.p2g.csv',
    params:
        thr_coaccess=config['methods']['celloracle']['thr_coaccess'],
        ext=config['methods']['celloracle']['ext']
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    shell:
        """
        set +e
        timeout $(({resources.runtime}-20))m bash -c \
        'Rscript workflow/scripts/mth/celloracle/p2g.R \
        {input.pre} \
        {input.csz} \
        {params.ext} \
        {output.pp} \
        {output.pc} && \
        python workflow/scripts/mth/celloracle/p2g.py \
        -d {input.pre} \
        -a {output.pp} \
        -c {output.pc} \
        -o {input.csz} \
        -t {params.thr_coaccess} \
        -p {output.out}'
        if [ $? -eq 124 ]; then
            awk 'BEGIN {{ print "cre,gene,score,pval" }}' > {output.out}
            touch {output.pp}
            touch {output.pc}
        fi
        """


rule tfb_celloracle:
    threads: 32
    singularity: 'workflow/envs/celloracle.sif'
    input:
        pre=lambda wildcards: map_rules('pre', wildcards.pre),
        p2g=lambda wildcards: map_rules('p2g', wildcards.p2g),
        gnm=rules.gen_genome_celloracle.output
    output: out='dts/{dat}/cases/{case}/runs/{pre}.{p2g}.celloracle.tfb.csv'
    params:
        fpr=config['methods']['celloracle']['fpr'],
        blen=config['methods']['celloracle']['blen'],
        tfb_thr=config['methods']['celloracle']['tfb_thr']
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    shell:
        """
        set +e
        timeout $(({resources.runtime}-20))m \
        python workflow/scripts/mth/celloracle/tfb.py \
        -d {input.pre} \
        -p {input.p2g} \
        -g {input.gnm} \
        -f {params.fpr} \
        -b {params.blen} \
        -t {params.tfb_thr} \
        -o {output.out}
        if [ $? -eq 124 ]; then
            awk 'BEGIN {{ print "cre,tf,score" }}' > {output.out}
        fi
        """


rule mdl_celloracle:
    threads: 32
    singularity: 'workflow/envs/celloracle.sif'
    input:
        pre=lambda wildcards: map_rules('pre', wildcards.pre),
        p2g=lambda wildcards: map_rules('p2g', wildcards.p2g),
        tfb=lambda wildcards: map_rules('tfb', wildcards.tfb),
    output:
        out='dts/{dat}/cases/{case}/runs/{pre}.{p2g}.{tfb}.celloracle.mdl.csv'
    params:
        a=config['methods']['celloracle']['a'],
        p=config['methods']['celloracle']['p'],
        n=config['methods']['celloracle']['n'],
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    shell:
        """
        export MKL_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        export NUMEXPR_NUM_THREADS=1
        set +e
        timeout $(({resources.runtime}-20))m \
        python workflow/scripts/mth/celloracle/mdl.py \
        -m {input.pre} \
        -g {input.p2g} \
        -t {input.tfb} \
        -a {params.a} \
        -p {params.p} \
        -n {params.n} \
        -o {output.out}
        if [ $? -eq 124 ]; then
            awk 'BEGIN {{ print "source,target,score,pval" }}' > {output.out}
        fi
        """


rule mdl_o_celloracle:
    threads: 32
    singularity: 'workflow/envs/celloracle.sif'
    input:
        mdata=rules.extract_case.output.mdata,
        csz=rules.gen_genome_celloracle.output,
    output:
        pp=temp(local('dts/{dat}/cases/{case}/runs/celloracle.src.peaks.csv')),
        pc=temp(local('dts/{dat}/cases/{case}/runs/celloracle.src.conns.csv')),
        out='dts/{dat}/cases/{case}/runs/o_celloracle.o_celloracle.o_celloracle.o_celloracle.mdl.csv',
    params:
        organism=lambda w: config['dts'][w.dat]['organism'],
        k=config['methods']['celloracle']['k'],
        thr_coaccess=config['methods']['celloracle']['thr_coaccess'],
        ext=config['methods']['celloracle']['ext'],
        fpr=config['methods']['celloracle']['fpr'],
        blen=config['methods']['celloracle']['blen'],
        tfb_thr=config['methods']['celloracle']['tfb_thr'],
        a=config['methods']['celloracle']['a'],
        p=config['methods']['celloracle']['p'],
        n=config['methods']['celloracle']['n'],
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'] * 2,
    shell:
        """
        export MKL_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        export NUMEXPR_NUM_THREADS=1
        set +e
        timeout $(({resources.runtime}-20))m bash -c \
        'Rscript workflow/scripts/mth/celloracle/src.R \
        {input.mdata} \
        {input.csz} \
        {params.ext} \
        {output.pp} \
        {output.pc} && \
        python workflow/scripts/mth/celloracle/src.py \
        -a {input.mdata} \
        -b {output.pp} \
        -c {output.pc} \
        -d {params.organism} \
        -e {params.thr_coaccess} \
        -f {params.fpr} \
        -g {params.blen} \
        -i {params.tfb_thr} \
        -j {params.a} \
        -k {params.p} \
        -l {params.n} \
        -m {params.k} \
        -n {output.out}'
        if [ $? -eq 124 ]; then
            awk 'BEGIN {{ print "source,target,score,pval" }}' > {output.out}
            touch {output.pp}
            touch {output.pc}
        fi
        """
