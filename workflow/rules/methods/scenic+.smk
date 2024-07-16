rule download_motifs:
    input:
        url_h = "https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
        url_h = "https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl"
    output:
        h = "aertslab/motifs-v10nr_clust/nr.hgnc-m0.001-o0.0.tbl",
        m = "aertslab/motifs-v10nr_clust/nr.mgi-m0.001-o0.0.tbl"
    shell:
        """
        wget -O {output.h} {input.url_h}cd en
        wget -O {output.m} {input.url_m}
        """
        


rule pre_scenicplus:
    input:
        cisTopic_obj_fname=config["input_data"]["cisTopic_obj_fname"],
        GEX_anndata_fname=config["input_data"]["GEX_anndata_fname"]
    output:
        config["output_data"]["combined_GEX_ACC_mudata"]
    params:
        bc_transform_func=lambda wildcards: config["params_data_preparation"]["bc_transform_func"]
    shell:
        """
        scenicplus prepare_data prepare_GEX_ACC \
            --cisTopic_obj_fname {input.cisTopic_obj_fname} \
            --GEX_anndata_fname {input.GEX_anndata_fname} \
            --out_file {output} \
            --bc_transform_func {params.bc_transform_func}
        """


motif_enrichment_cistarget
download_genome_annotations
motif_enrichment_dem, prepare_menr, get_search_space, region_to_gene  
tf_to_gene, eGRN_direct, eGRN_extended