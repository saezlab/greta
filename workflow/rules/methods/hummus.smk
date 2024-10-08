localrules: pre_hummus


rule pre_hummus:
    input:
        mdata='datasets/{dataset}/cases/{case}/mdata.h5mu',
    output:
        out='datasets/{dataset}/cases/{case}/runs/hummus.pre.h5mu'
    params:
        organism=lambda w: config['datasets'][w.dataset]['organism'],
    shell:
        """
        cp {input.mdata} {output.out}
        """


rule prior_hummus:
    threads: 16
    input:
        pre=lambda wildcards: map_rules('pre', wildcards.pre),
        h=rules.download_granges.output.h,
        m=rules.download_granges.output.m,
    output:
        tf_list=temp(local('datasets/{dataset}/cases/{case}/runs/{pre}.tfs.list.hummus.tsv')),
        rna_network=temp(local('datasets/{dataset}/cases/{case}/runs/{pre}.grnboost2.csv')),
        atac_network=temp(local('datasets/{dataset}/cases/{case}/runs/{pre}.atacnet.csv')),
        hummus_object=temp(local('datasets/{dataset}/cases/{case}/runs/{pre}.hummus_object.RDS')),
    params:
        organism=lambda w: config['datasets'][w.dataset]['organism'],
        grn_number_edges=50000,
        tf_layer_method=None,
        n_cores=16
    singularity:
        'workflow/envs/hummus.sif'
    resources:
        mem_mb=128000,
    shell:
        """
        export OMP_NUM_THREADS={params.n_cores}
        export MKL_NUM_THREADS={params.n_cores}
        export NUMEXPR_NUM_THREADS={params.n_cores}

        Rscript workflow/scripts/methods/hummus/prior_hummus_tf_infos.R \
        {input.pre} \
        {params.organism} \
        {output.tf_list}

        python workflow/scripts/methods/hummus/prior_hummus.py \
        -d {input.pre} \
        -r {output.rna_network} \
        -a {output.atac_network} \
        -c {params.n_cores} \
        -n {output.tf_list} \
        -o {params.organism}
        echo "Nets and circe done"
        echo "Starting last step"
        Rscript workflow/scripts/methods/hummus/prior_hummus.R \
        {input.pre} \
        {params.organism} \
        {input.h} \
        {input.m} \
        {output.rna_network} \
        {output.atac_network} \
        {output.hummus_object} \
        {params.grn_number_edges} \
        {params.tf_layer_method}
        """


rule p2g_hummus:
    threads: 32
    input:
        hummus_object=rules.prior_hummus.output.hummus_object
    singularity:
        'workflow/envs/hummus.sif'
    output:
        out='datasets/{dataset}/cases/{case}/runs/{pre}.hummus.p2g.csv',
        multilayer_f=temp(directory('datasets/{dataset}/cases/{case}/runs/{pre}.multilayer')),
    params:
        organism=lambda w: config['datasets'][w.dataset]['organism'],
        ext=900000,
        n_cores=32
    shell:
        """
        Rscript workflow/scripts/methods/hummus/p2g.R \
        {input.hummus_object} \
        {params.organism} \
        {params.ext} \
        {params.n_cores} \
        {output.out} \
        {output.multilayer_f} \
        """


rule tfb_hummus:
    threads: 32
    input:
        p2g=lambda wildcards: map_rules('p2g', wildcards.p2g),
        hummus_object=rules.prior_hummus.output.hummus_object,
    singularity:
        'workflow/envs/hummus.sif'
    output:
        out='datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.hummus.tfb.csv',
        multilayer_f=temp(directory('datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.multilayer')),
    params:
        organism=lambda w: config['datasets'][w.dataset]['organism'],
        n_cores=32,
        p2g=lambda w: w.p2g
    shell:
        """
        Rscript workflow/scripts/methods/hummus/tfb.R \
        {input.hummus_object} \
        {input.p2g} \
        {params.organism} \
        {params.n_cores} \
        {output.out} \
        {output.multilayer_f} \
        {params.p2g}
        """


rule mdl_hummus:
    threads: 32
    input:
        p2g=lambda wildcards: map_rules('p2g', wildcards.p2g),
        tfb=lambda wildcards: map_rules('tfb', wildcards.tfb),
        hummus_object=rules.prior_hummus.output.hummus_object,
    singularity:
        'workflow/envs/hummus.sif'
    output:
        out='datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.hummus.mdl.csv',
        multilayer_f=temp(directory('datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.multilayer')),
    params:
        n_cores=32,
        organism=lambda w: config['datasets'][w.dataset]['organism'],
        p2g=lambda w: w.p2g,
        tfb=lambda w: w.tfb,
    shell:
        """
        Rscript workflow/scripts/methods/hummus/mdl.R \
        {input.hummus_object} \
        {input.p2g} \
        {input.tfb} \
        {params.n_cores} \
        {output.out} \
        {output.multilayer_f} \
        {params.p2g} \
        {params.tfb} \
        {params.organism}
        """
