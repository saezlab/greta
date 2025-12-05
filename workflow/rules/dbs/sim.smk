localrules: download_sim

rule download_sim:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        img='workflow/envs/gretabench.sif',
    output: directory('dbs/sim')
    params:
        url=config['url_sim'],
    shell:
        """
        path_sim=$(dirname {output})
        path_tmp=$path_sim/tmp.zip
        wget --no-verbose "{params.url}" -O $path_tmp
        python -c "import zipfile, os; os.makedirs('$path_sim', exist_ok=True); [zipfile.ZipFile('$path_tmp').extract(m, '$path_sim') for m in zipfile.ZipFile('$path_tmp').namelist()]"
        mv $path_sim/simulations {output}
        rm -rf $path_tmp
        for f in {output}/seed_*; do
            base=$(basename "$f")
            path_dataset={output}/$base
            python workflow/scripts/dbs/sim/pre.py $path_dataset
            python workflow/scripts/dts/extract_case.py \
            -i $path_dataset/mdata.h5mu -c all -s 0 -d 0 -g 20000 -r 80000 -t None -o $path_dataset/mdata.h5mu
        done
        """


rule sim_pearson:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: rules.download_sim.output
    output: 'dts/sim/seed_{sim_num}/pearson.csv'
    params:
        mode=config['methods']['pearson']['mode'],
        thr_r2=config['methods']['pearson']['thr_r2'],
    shell:
        """
        python workflow/scripts/mth/correlation.py \
        -i {input}/seed_{wildcards.sim_num}/mdata.h5mu \
        -t {input}/seed_{wildcards.sim_num}/tfs.txt \
        -m {params.mode} \
        -r {params.thr_r2} \
        -o {output}
        """


rule sim_spearman:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: rules.download_sim.output
    output: 'dts/sim/seed_{sim_num}/spearman.csv'
    params:
        mode=config['methods']['spearman']['mode'],
        thr_r2=config['methods']['spearman']['thr_r2'],
    shell:
        """
        python workflow/scripts/mth/correlation.py \
        -i {input}/seed_{wildcards.sim_num}/mdata.h5mu \
        -t {input}/seed_{wildcards.sim_num}/tfs.txt \
        -m {params.mode} \
        -r {params.thr_r2} \
        -o {output}
        """


rule sim_random:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: rules.download_sim.output
    output: 'dts/sim/seed_{sim_num}/random.csv'
    params:
        g_perc=config['methods']['random']['g_perc'],
        scale=config['methods']['random']['scale'],
        tf_g_ratio=config['methods']['random']['tf_g_ratio'],
        w_size=config['methods']['random']['w_size'] // 2,
        seed=config['methods']['random']['seed'],
    shell:
        """
        python workflow/scripts/mth/random/grn.py \
        -i {input}/seed_{wildcards.sim_num}/mdata.h5mu \
        -t {input}/seed_{wildcards.sim_num}/tfs.txt \
        -g {params.g_perc} \
        -n {params.scale} \
        -r {params.tf_g_ratio} \
        -w {params.w_size} \
        -s {params.seed} \
        -o {output}
        """


rule sim_grnboost:
    threads: 2
    singularity: 'workflow/envs/scenicplus.sif'
    input:
        img='workflow/envs/scenicplus.sif',
        inp=rules.download_sim.output
    output:
        out='dts/sim/seed_{sim_num}/grnboost.csv'
    shell:
        """
        # Create Loom file
        path_tmp=$(dirname {output.out})
        path_tmp=$path_tmp/tmp_grnboost
        path_gex=$path_tmp/gex.loom
        path_adj=$path_tmp/adj.tsv
        rm -rf $path_tmp
        mkdir -p $path_tmp
        python workflow/scripts/mth/scenic/loom.py \
        -i {input.inp}/seed_{wildcards.sim_num}/mdata.h5mu \
        -o $path_gex
        echo "Created loom"
        # Run grnboost2 GRN
        arboreto_with_multiprocessing.py $path_gex {input.inp}/seed_{wildcards.sim_num}/tfs.txt -o $path_adj --num_workers {threads} --seed 42
        echo "Generated adj"
        # Process GRN
        python workflow/scripts/mth/scenic/process_grn.py \
        -g $path_adj \
        -o {output.out}
        rm -rf $path_tmp
        """


rule sim_celloracle:
    threads: 1
    singularity: 'workflow/envs/celloracle.sif'
    input:
        img='workflow/envs/celloracle.sif',
        inp=rules.download_sim.output
    output:
        out='dts/sim/seed_{sim_num}/celloracle.csv'
    params:
        a=config['methods']['celloracle']['a'],
        p=config['methods']['celloracle']['p'],
        n=config['methods']['celloracle']['n'],
    shell:
        """
        export MKL_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        export NUMEXPR_NUM_THREADS=1
        python workflow/scripts/mth/celloracle/mdl.py \
        -m {input.inp}/seed_{wildcards.sim_num}/mdata.h5mu \
        -g {input.inp}/seed_{wildcards.sim_num}/p2g.csv \
        -t {input.inp}/seed_{wildcards.sim_num}/tfb.csv \
        -a {params.a} \
        -p {params.p} \
        -n {params.n} \
        -o {output}
        """


rule sim_pando:
    threads: 1
    singularity: 'workflow/envs/pando.sif'
    input:
        img='workflow/envs/pando.sif',
        inp=rules.download_sim.output,
    output:
        out='dts/sim/seed_{sim_num}/pando.csv'
    params:
        thr_corr=0,
        p_thresh=0.1,
        rsq_thresh=0,
        nvar_thresh=0,
        min_genes_per_module=config['methods']['pando']['min_genes_per_module'],
    shell:
        """
        Rscript workflow/scripts/mth/pando/mdl.R \
        {input.inp}/seed_{wildcards.sim_num}/mdata.h5mu \
        {input.inp}/seed_{wildcards.sim_num}/p2g.csv \
        {input.inp}/seed_{wildcards.sim_num}/tfb.csv \
        {params.thr_corr} \
        {params.p_thresh} \
        {params.rsq_thresh} \
        {params.nvar_thresh} \
        {params.min_genes_per_module} \
        {threads} \
        {output.out}
        """

rule sim_figr:
    threads: 1
    singularity: 'workflow/envs/figr.sif'
    input:
        img='workflow/envs/figr.sif',
        inp=rules.download_sim.output,
    output:
        out='dts/sim/seed_{sim_num}/figr.csv'
    params:
        cellK=config['methods']['figr']['cellK'],
        thr_score=0.5,
    shell:
        """
        Rscript workflow/scripts/mth/figr/mdl.R \
        {input.inp}/seed_{wildcards.sim_num}/mdata.h5mu \
        {input.inp}/seed_{wildcards.sim_num}/p2g.csv \
        {input.inp}/seed_{wildcards.sim_num}/tfb.csv \
        {params.cellK} \
        {params.thr_score} \
        {threads} \
        {output.out}
        """
