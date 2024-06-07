# List of rules
# 1. pre_hummus
# 2. peak to genes
    # 2.1 tf_motifs and tf_list (use internet)
    # 2.2 grnboost2
    # 2.3 create cicero metacells
    # 2.4 create cicero network
    # 2.5 create hummus object and add computed networks
    # 2.6 p2g_hummus

# 3. tf binding
    # 3.1 tfb_hummus (update p2g bipartite in hummus object and compute tfb_hummus)

rule pre_hummus:
    input:
        d='datasets/{dataset}/cases/{case}/mdata.h5mu',
        #h='gdata/grandes/EnsDb.Hsapiens.v86.RDS',
    singularity:
        'workflow/envs/gretabench.sif'
    benchmark:
        'benchmarks/{dataset}.{case}.hummus.pre.txt'
    output:
        d='datasets/{dataset}/cases/{case}/runs/hummus.pre.h5mu'
    params:
        organism=lambda w: config['datasets'][w.dataset]['organism'],
    shell:
        """
        python workflow/scripts/methods/hummus/pre.py \
        -i {input.d} \
        -o {output.d}
        """

rule prior_hummus:
    input:
        d='datasets/{dataset}/cases/{case}/runs/{pre}.pre.h5mu',
        h = 'gdata/granges/hg38_ensdb_v86.csv',
        m = 'gdata/granges/mm10_ensdb_v79.csv',
    output:
        tf_list = temp('datasets/{dataset}/cases/{case}/runs/{pre}.tfs.list.hummus.tsv'),
        rna_network = temp('datasets/{dataset}/cases/{case}/runs/{pre}.grnboost2.csv'),
        atac_network = temp('datasets/{dataset}/cases/{case}/runs/{pre}.atacnet.csv'),
        hummus_object = temp('datasets/{dataset}/cases/{case}/runs/{pre}.hummus_object.RDS'),
    params:
        organism=lambda w: config['datasets'][w.dataset]['organism'],
        distance_threshold_atacnet = 500000,
        n_cores = 32,
    singularity:
        'workflow/envs/hummus.sif'
    shell:
        """
        Rscript workflow/scripts/methods/hummus/get_tf_infos.R \
        {input.d} \
        {params.organism} \
        {output.tf_list}

        python workflow/scripts/methods/hummus/prior_hummus.py \
        -d {input.d} \
        -r {output.rna_network} \
        -a {output.atac_network} \
        -t {params.distance_threshold_atacnet} \
        -c {params.n_cores} \
        -n {output.tf_list}

        Rscript workflow/scripts/methods/hummus/get_tf_infos.R \
        {params.organism} \
        {input.h} \
        {input.m} \
        {output.rna_network} \
        {output.atac_network} \
        {output.hummus_object}
        """


#rule grnboost2:
#    input:
#        d='datasets/{dataset}/cases/{case}/runs/{pre}.pre.h5mu',
#        tf_list='datasets/{dataset}/cases/{case}/runs/{pre}.tfs.list.hummus.tsv'
#    singularity:
#        'workflow/envs/grnboost2.sif'
#    benchmark:
#        'benchmarks/{dataset}.{case}.{pre}.hummus.grnboost2.txt'
#    params:
#        n_cores = 30
#    output:
#        rna_loom = temp('datasets/{dataset}/cases/{case}/runs/{pre}.rna.loom'),
#        rna_network = temp(local('datasets/{dataset}/cases/{case}/runs/{pre}.grnboost2.csv'))
#    shell:
#        """
#        #add prepare rna as loom TODO add/change to d instead of rna
#        python workflow/scripts/methods/hummus/make_loom_from_h5mu.py \
#        -f {input.d} \
#        -o {output.rna_loom}
#
#        #add compute grn_boost_network from loom
#        pyscenic grn \
#        {output.rna_loom} \
#        {input.tf_list} \
#        -o {output.rna_network} \
#        --num_workers {params.n_cores} \
#        --seed 100 --sparse
#        """

#rule save_pseudocells_cicero:
#    input:
#        atac_path = "datasets/{dataset}/cases/{case}/runs/{pre}.pre.h5mu"
#    output:
#        cicero_cds_tsv = "datasets/{dataset}/cases/{case}/runs/{pre}.cicero_cds.tsv"
#    params:
#        number_cells_per_clusters = 50,
#    singularity:
#        "workflow/envs/hummus.sif"
#    script:
#        "../../../workflow/scripts/methods/hummus/cicero_save_cds.R"

#rule run_atac_networks:
#    input:
##        cicero_cds = "datasets/{dataset}/cases/{case}/runs/{pre}.cicero_cds.tsv",
#        mudata = "datasets/{dataset}/cases/{case}/runs/{pre}.pre.h5mu"
#    output:
#        atacnet_network = "datasets/{dataset}/cases/{case}/runs/{pre}.atacnet.csv"
#    params:
#        distance_threshold = 500000,
##    benchmark:
# #       "benchmarks/{dataset}.{case}.{pre}.atacnet_network.txt"
#    singularity:
#        "workflow/envs/hummus.sif"
#    script:
#        "../../../workflow/scripts/methods/hummus/atac_network.py"
#
#rule initiate_multilayer:
#    input:
#        d='datasets/{dataset}/cases/{case}/runs/{pre}.pre.h5mu',
#        h='gdata/granges/hg38_ensdb_v86.csv',
#        m='gdata/granges/mm10_ensdb_v79.csv',
#        rna_network = local('datasets/{dataset}/cases/{case}/runs/{pre}.grnboost2.csv'),
#        atac_network = local('datasets/{dataset}/cases/{case}/runs/{pre}.atacnet.csv')
#    singularity:
#        'workflow/envs/hummus.sif'
#    benchmark:
#        'benchmarks/{dataset}.{case}.{pre}.hummus.p2g.txt'
#    output:
#        # multilayer_f indictaed through multiplex subfolder that won't change
#        # bipartites not given since they can be change with other p2g and tfb outputs
#        hummus_object = "datasets/{dataset}/cases/{case}/runs/{pre}.hummus_object.RDS",
#    params:
#        organism=lambda w: config['datasets'][w.dataset]['organism'],
#        n_cores = 32
#    shell:
#        """
#        # add Integrate networks into multilayer object, and connect them
#        Rscript workflow/scripts/methods/hummus/initiate_multilayer.R \
#        {input.d} \
#        {params.organism} \
#        {input.h} \
#        {input.m} \
#        {input.rna_network} \
#        {input.atac_network} \
#        {output.hummus_object} \
#        """
#

rule p2g_hummus:
    input:
        hummus_object = "datasets/{dataset}/cases/{case}/runs/{pre}.hummus_object.RDS"
    singularity:
        'workflow/envs/hummus.sif'
    benchmark:
        'benchmarks/{dataset}.{case}.{pre}.hummus.p2g.txt'
    output:
        p2g = 'datasets/{dataset}/cases/{case}/runs/{pre}.hummus.p2g.csv',
#        loom_rna = temp('datasets/{dataset}/cases/{case}/runs/{pre}.hummus.rna.loom'),
#        atac_network = temp(local('datasets/{dataset}/cases/{case}/runs/{pre}.hummus.atacnet.csv')),
        # multilayer_f indictaed through multiplex subfolder that won't change
        # bipartites not given since they can be change with other p2g and tfb outputs
        multilayer_f = temp(directory("datasets/{dataset}/cases/{case}/runs/{pre}.multilayer")),
    params:
        organism=lambda w: config['datasets'][w.dataset]['organism'],
        ext=900000,
        n_cores = 32
    shell:
        """
        # add Integrate networks into multilayer object, and connect them
        Rscript workflow/scripts/methods/hummus/p2g.R \
        {input.hummus_object} \
        {params.organism} \
        {params.ext} \
        {params.n_cores} \
        {output.p2g} \
        {output.multilayer_f} \
        """


rule tfb_hummus:
    input:
        #d='datasets/{dataset}/cases/{case}/runs/{pre}.pre.h5mu',
        p2g = 'datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.p2g.csv',
        hummus_object = 'datasets/{dataset}/cases/{case}/runs/{pre}.hummus_object.RDS'
    singularity:
        'workflow/envs/hummus.sif'
    benchmark:
        'benchmarks/{dataset}.{case}.{pre}.{p2g}.hummus.tfb.txt'
    output:
        tfb = 'datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.hummus.tfb.csv',
        multilayer_f = temp(directory("datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.multilayer")),
    params:
        organism=lambda w: config['datasets'][w.dataset]['organism'],
        n_cores = 32
    shell:
        """
        Rscript workflow/scripts/methods/hummus/tfb.R \
        {input.hummus_object} \
        {input.p2g} \
        {params.organism} \
        {params.n_cores} \
        {output.tfb} \
        {output.multilayer_f} \
        """

rule mdl_hummus:
    input:
        #d='datasets/{dataset}/cases/{case}/runs/{pre}.pre.h5mu',
        p2g = 'datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.p2g.csv',
        tfb = 'datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.tfb.csv',
        hummus_object = 'datasets/{dataset}/cases/{case}/runs/{pre}.hummus_object.RDS'
    singularity:
        'workflow/envs/hummus.sif'
    benchmark:
        'benchmarks/{dataset}.{case}.{pre}.{p2g}.{tfb}.hummus.mdl.txt'
    output:
        mdl = 'datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.hummus.mdl.csv',
        multilayer_f = temp(directory("datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.multilayer")),
    params:
        n_cores = 32
    shell:
        """
        Rscript workflow/scripts/methods/hummus/mdl.R \
        {input.hummus_object} \
        {input.p2g} \
        {input.tfb} \
        {params.n_cores} \
        {output.mdl} \
        {output.multilayer_f} \
        """
