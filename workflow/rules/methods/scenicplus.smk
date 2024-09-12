rule download_motifs:
    params:
        url_h = "https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl",
        url_m = "https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl"
    output:
        h = "aertslab/motifs-v10nr_clust/nr.hgnc-m0.001-o0.0.tbl",
        m = "aertslab/motifs-v10nr_clust/nr.mgi-m0.001-o0.0.tbl"
    shell:
        """
        wget -O {output.h} {params.url_h}
        wget -O {output.m} {params.url_m}
        """

rule download_gene_annotations:
    output:
        h = "aertslab/genomes/hg38/hg38_ensdb_v86.csv",
        m = "aertslab/genomes/mm10/mm10_ensdb_v79.csv"
    singularity:
        "workflow/envs/scenicplus.sif"
    shell:
        """
        python workflow/scripts/methods/scenicplus/download_gene_annot.py \
        -j {output.h} \
        -m {output.m} 
        """

rule download_chrom_sizes:
    output:
        m = "aertslab/genomes/mm10/mm10.chrom.sizes",
        h = "aertslab/genomes/hg38/hg38.chrom.sizes"
    shell:
        """
        wget -O {output.m} http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
        wget -O {output.h} http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes
        """


rule download_cistarget:
    output:
        human_rankings = "aertslab/cistarget/human_motif_SCREEN.regions_vs_motifs.rankings.feather",
        human_scores = "aertslab/cistarget/human_motif_SCREEN.regions_vs_motifs.scores.feather",
        mouse_rankings = "aertslab/cistarget/mouse_motif_SCREEN.regions_vs_motifs.rankings.feather",
        mouse_scores = "aertslab/cistarget/mouse_motif_SCREEN.regions_vs_motifs.scores.feather"
    params:
        human_rankings_url = "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/screen/mc_v10_clust/region_based/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather",
        human_scores_url = "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/screen/mc_v10_clust/region_based/hg38_screen_v10_clust.regions_vs_motifs.scores.feather",
        mouse_rankings_url = "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/screen/mc_v10_clust/region_based/mm10_screen_v10_clust.regions_vs_motifs.rankings.feather",
        mouse_scores_url = "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/screen/mc_v10_clust/region_based/mm10_screen_v10_clust.regions_vs_motifs.scores.feather"
    shell:
        """
        wget -O {output.human_rankings} {params.human_rankings_url}
        wget -O {output.human_scores} {params.human_scores_url}
        wget -O {output.mouse_rankings} {params.mouse_rankings_url}
        wget -O {output.mouse_scores} {params.mouse_scores_url}
        """


rule pre_scenicplus:
    input:
        frags = 'datasets/{dataset}/smpl.frags.tsv.gz',
        mudata = 'datasets/{dataset}/cases/{case}/mdata.h5mu',
        chrom_sizes_m = "aertslab/genomes/mm10/mm10.chrom.sizes",
        chrom_sizes_h = "aertslab/genomes/hg38/hg38.chrom.sizes",
        annot_m = "aertslab/genomes/mm10/mm10_ensdb_v79.csv",
        annot_h = "aertslab/genomes/hg38/hg38_ensdb_v86.csv"
    output:
        d = 'datasets/{dataset}/cases/{case}/runs/scenicplus.pre.h5mu',
    singularity:
        'workflow/envs/scenicplus.sif'
    params:
        tmp_scenicplus = temp(directory(local('datasets/{dataset}/cases/{case}/runs/scenicplus_tmp'))), 
        organism=lambda w: config['datasets'][w.dataset]['organism'],
        n_cores = 32
    shell:
        """
        python workflow/scripts/methods/scenicplus/pre.py \
        -f {input.frags} \
        -i {input.mudata} \
        -m {input.chrom_sizes_m} \
        -j {input.chrom_sizes_h} \
        -o {output.d} \
        -t {params.tmp_scenicplus} \
        -g {params.organism} \
        -n {params.n_cores} \
        -a {input.annot_m} \
        -b {input.annot_h}
        """

rule p2g_scenicplus:
    input:
        data='datasets/{dataset}/cases/{case}/runs/{pre}.pre.h5mu',
        chrom_sizes_m = "aertslab/genomes/mm10/mm10.chrom.sizes",
        chrom_sizes_h = "aertslab/genomes/hg38/hg38.chrom.sizes"
    singularity:
        'workflow/envs/scenicplus.sif'
    params:
        temp_dir = temp(directory(local('datasets/{dataset}/cases/{case}/runs/scenicplus_tmp/{pre}'))),
        n_cores = 32,
        organism=lambda w: config['datasets'][w.dataset]['organism'],
        # Minimum and maximum (up until another gene) number of bps upstream to include in the search space.
        # cf : scenicplus.data_wrangling.gene_search_space.get_search_space
        search_space_upstream = "1000 150000",
        search_space_downstream = "1000 150000",
        search_space_extend_tss = "10 10",
        remove_promoters = False,  # Fixed to False in the snakemake pipeline
        use_gene_boundaries = False,  # Fixed to False in the snakemake pipeline
        region_to_gene_importance_method = "GBM",
        region_to_gene_correlation_method = "SR",
    benchmark:
        'benchmarks/{dataset}.{case}.{pre}.scenicplus.p2g.txt'
    output:
        pg='datasets/{dataset}/cases/{case}/runs/{pre}.scenicplus.p2g.csv',
    shell:
        """
        python workflow/scripts/methods/scenicplus/p2g.py \
        -i {input.data} \
        -c {params.n_cores} \
        -t {params.temp_dir} \
        -m {input.chrom_sizes_m} \
        -j {input.chrom_sizes_h} \
        -o {output.pg} \
        -u "{params.search_space_upstream}" \
        -d "{params.search_space_downstream}" \
        -e "{params.search_space_extend_tss}" \
        -r {params.remove_promoters} \
        -g {params.organism} \
        -z {params.use_gene_boundaries} \
        -p {params.region_to_gene_importance_method} \
        -s {params.region_to_gene_correlation_method}
        """

rule tfb_scenicplus:
    input:
        raw = "datasets/{dataset}/cases/{case}/mdata.h5mu",
        pre = "datasets/{dataset}/cases/{case}/runs/{pre}.pre.h5mu",
        p2g = "datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.p2g.csv",
        cistarget_rankings_human = "aertslab/cistarget/human_motif_SCREEN.regions_vs_motifs.rankings.feather",
        cistarget_scores_human = "aertslab/cistarget/human_motif_SCREEN.regions_vs_motifs.scores.feather",
        path_to_motif_annotations_human = "aertslab/motifs-v10nr_clust/nr.hgnc-m0.001-o0.0.tbl",
        path_to_motif_annotations_mouse = "aertslab/motifs-v10nr_clust/nr.mgi-m0.001-o0.0.tbl"
    params:
        n_cores = 32,
        organism=lambda w: config['datasets'][w.dataset]['organism'],
    output:
        cistarget_results = temp("datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.scenicplus.cistarget.hdf5"),
        dem_results = temp("datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.scenicplus.dem.hdf5"),
        annotation_direct_path = temp("datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.scenicplus.annotation_direct.h5ad"),
        annotation_extended_path = temp("datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.scenicplus.annotation_extended.h5ad"),
        tf_names_path = temp("datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.scenicplus.tf_names.txt"),
        tfb = "datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.scenicplus.tfb.csv"
    singularity:
        'workflow/envs/scenicplus.sif'
    shell:
        """
        python workflow/scripts/methods/scenicplus/tfb.py \
        -i {input.pre} \
        -p {input.p2g} \
        -d {input.raw} \
        -r {input.cistarget_rankings_human} \
        -s {input.cistarget_scores_human} \
        -g {params.organism} \
        -t {output.cistarget_results} \
        -u {output.dem_results} \
        -o {output.tfb} \
        -c {params.n_cores} \
        --annotation_direct_path {output.annotation_direct_path}\
        --annotation_extended_path {output.annotation_extended_path}\
        --tf_names_path {output.tf_names_path}\
        --path_to_motif_annotations_human {input.path_to_motif_annotations_human}\
        --path_to_motif_annotations_mouse {input.path_to_motif_annotations_mouse}
        """
# motif_enrichment_cistarget
# download_genome_annotations
# motif_enrichment_dem, prepare_menr, get_search_space, region_to_gene  
# tf_to_gene, eGRN_direct, eGRN_extended

rule mdl_scenicplus:
    input:
        mdata = "datasets/{dataset}/cases/{case}/runs/{pre}.pre.h5mu",
        p2g = "datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.p2g.csv",
        tfb = "datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.tfb.csv",
        cistarget_rankings_human = "aertslab/cistarget/human_motif_SCREEN.regions_vs_motifs.rankings.feather",
        cistarget_rankings_mouse = "aertslab/cistarget/mouse_motif_SCREEN.regions_vs_motifs.rankings.feather",
        #tf_names_path = "datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.scenicplus.tf_names.txt",
    output:
#        tf_gene_prior_output = temp("datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.scenicplus.tf_gene_prior.csv"),
#        eRegulon_output = temp("datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.scenicplus.eRegulon.csv"),
        mdl = "datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.scenicplus.mdl.csv"
    params:
        method_mdl = "GBM",
        n_cores = 32,
        organism=lambda w: config['datasets'][w.dataset]['organism'],
        temp_dir = "/tmp",
        order_regions_to_genes_by="importance",
        order_TFs_to_genes_by="importance",
        gsea_n_perm=1000,
        quantile_thresholds_region_to_gene="0.85 0.90 0.95",
        top_n_regionTogenes_per_gene="5 10 15",
        top_n_regionTogenes_per_region="",
        min_regions_per_gene=0,
        rho_threshold=0.05,
        min_target_genes=10,
    singularity:
        'workflow/envs/scenicplus.sif'
    shell:
        """
        python workflow/scripts/methods/scenicplus/mdl.py \
        -i {input.mdata} \
        -p {input.p2g} \
        -t {input.tfb} \
        -l {input.cistarget_rankings_human} \
        -k {input.cistarget_rankings_mouse} \
        -o {output.mdl} \
        -m {params.method_mdl} \
        -c {params.n_cores} \
        -g {params.organism} \
        --temp_dir {params.temp_dir} \
        --order_regions_to_genes_by {params.order_regions_to_genes_by} \
        --order_TFs_to_genes_by {params.order_TFs_to_genes_by} \
        --gsea_n_perm {params.gsea_n_perm} \
        --quantile_thresholds_region_to_gene {params.quantile_thresholds_region_to_gene} \
        --top_n_regionTogenes_per_gene {params.top_n_regionTogenes_per_gene} \
        --top_n_regionTogenes_per_region {params.top_n_regionTogenes_per_region} \
        --min_regions_per_gene {params.min_regions_per_gene} \
        --rho_threshold {params.rho_threshold} \
        --min_target_genes {params.min_target_genes}
        """

