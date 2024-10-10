localrules: download_motifs, download_gene_annotations, download_cistarget


rule download_motifs:
    params:
        url_h="https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl",
        url_m="https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl"
    output:
        h="aertslab/motifs-v10nr_clust/nr.hgnc-m0.001-o0.0.tbl",
        m="aertslab/motifs-v10nr_clust/nr.mgi-m0.001-o0.0.tbl"
    shell:
        """
        wget -O {output.h} {params.url_h}
        wget -O {output.m} {params.url_m}
        """


rule download_gene_annotations:
    output:
        h="aertslab/genomes/hg38/hg38_ensdb_v86.csv",
        m="aertslab/genomes/mm10/mm10_ensdb_v79.csv"
    singularity:
        "workflow/envs/scenicplus.sif"
    shell:
        """
        python workflow/scripts/methods/scenicplus/download_gene_annot.py \
        -j {output.h} \
        -m {output.m} 
        """


rule download_cistarget:
    output:
        human_rankings="aertslab/cistarget/human_motif_SCREEN.regions_vs_motifs.rankings.feather",
        human_scores="aertslab/cistarget/human_motif_SCREEN.regions_vs_motifs.scores.feather",
        mouse_rankings="aertslab/cistarget/mouse_motif_SCREEN.regions_vs_motifs.rankings.feather",
        mouse_scores="aertslab/cistarget/mouse_motif_SCREEN.regions_vs_motifs.scores.feather"
    params:
        human_rankings_url="https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/screen/mc_v10_clust/region_based/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather",
        human_scores_url="https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/screen/mc_v10_clust/region_based/hg38_screen_v10_clust.regions_vs_motifs.scores.feather",
        mouse_rankings_url="https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/screen/mc_v10_clust/region_based/mm10_screen_v10_clust.regions_vs_motifs.rankings.feather",
        mouse_scores_url="https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/screen/mc_v10_clust/region_based/mm10_screen_v10_clust.regions_vs_motifs.scores.feather"
    shell:
        """
        wget -O {output.human_rankings} {params.human_rankings_url}
        wget -O {output.human_scores} {params.human_scores_url}
        wget -O {output.mouse_rankings} {params.mouse_rankings_url}
        wget -O {output.mouse_scores} {params.mouse_scores_url}
        """


rule pre_scenicplus:
    threads: 32
    input:
        frags=list_frags_files,
        mdata=rules.extract_case.output.mdata,
        chrom_sizes_h=rules.download_genomesizes.output.hg38,
        chrom_sizes_m=rules.download_genomesizes.output.mm10,
        annot_m=rules.download_gene_annotations.output.m,
        annot_h=rules.download_gene_annotations.output.h,
    output:
        out='datasets/{dataset}/cases/{case}/runs/scenicplus.pre.h5mu',
    singularity:
        'workflow/envs/scenicplus.sif'
    params:
        tmp_scenicplus=temp(directory(local('datasets/{dataset}/cases/{case}/runs/scenicplus_tmp'))), 
        ray_tmp_dir="/tmp",
        organism=lambda w: config['datasets'][w.dataset]['organism'],
        n_cores=32
    shell:
        """
        python workflow/scripts/methods/scenicplus/pre.py \
        -f {input.frags} \
        -i {input.mdata} \
        -m {input.chrom_sizes_m} \
        -j {input.chrom_sizes_h} \
        -o {output.out} \
        -t {params.tmp_scenicplus} \
        -g {params.organism} \
        -n {params.n_cores} \
        -a {input.annot_m} \
        -b {input.annot_h} \
        --ray_tmp_dir '{params.ray_tmp_dir}'
        """


rule p2g_scenicplus:
    threads: 32
    input:
        pre=lambda wildcards: map_rules('pre', wildcards.pre),
        chrom_sizes_h=rules.download_genomesizes.output.hg38,
        chrom_sizes_m=rules.download_genomesizes.output.mm10,
    singularity:
        'workflow/envs/scenicplus.sif'
    params:
        temp_dir=temp(directory(local('datasets/{dataset}/cases/{case}/runs/scenicplus_tmp/{pre}'))),
        n_cores=32,
        organism=lambda w: config['datasets'][w.dataset]['organism'],
        ext=config['methods']['scenicplus']['ext'],
        min_dist=config['methods']['scenicplus']['min_dist'],
        search_space_extend_tss=config['methods']['scenicplus']['search_space_extend_tss'],
        remove_promoters=config['methods']['scenicplus']['remove_promoters'],
        use_gene_boundaries=config['methods']['scenicplus']['use_gene_boundaries'],
        region_to_gene_importance_method=config['methods']['scenicplus']['region_to_gene_importance_method'],
        region_to_gene_correlation_method=config['methods']['scenicplus']['region_to_gene_correlation_method'],
    output:
        out='datasets/{dataset}/cases/{case}/runs/{pre}.scenicplus.p2g.csv',
    shell:
        """
        python workflow/scripts/methods/scenicplus/p2g.py \
        -i {input.pre} \
        -c {params.n_cores} \
        -t {params.temp_dir} \
        -m {input.chrom_sizes_m} \
        -j {input.chrom_sizes_h} \
        -o {output.out} \
        -u "{params.ext}" \
        -d "{params.min_dist}" \
        -e "{params.search_space_extend_tss}" \
        -r {params.remove_promoters} \
        -g {params.organism} \
        -z {params.use_gene_boundaries} \
        -p {params.region_to_gene_importance_method} \
        -s {params.region_to_gene_correlation_method}
        """


rule tfb_scenicplus:
    threads: 32
    input:
        mdata=rules.extract_case.output.mdata,
        pre=lambda wildcards: map_rules('pre', wildcards.pre),
        p2g=lambda wildcards: map_rules('p2g', wildcards.p2g),
        cistarget_rankings_human=rules.download_cistarget.output.human_rankings,
        cistarget_scores_human=rules.download_cistarget.output.human_scores,
        path_to_motif_annotations_human=rules.download_motifs.output.h,
        path_to_motif_annotations_mouse=rules.download_motifs.output.m,
    params:
        n_cores=32,
        organism=lambda w: config['datasets'][w.dataset]['organism'],
        temp_dir=temp(directory(local("datasets/{dataset}/cases/{case}/runs/{pre}.scenicplus_tmp"))),
        ray_tmp_dir="/tmp",
    output:
        cistarget_results=temp("datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.scenicplus.cistarget.hdf5"),
        dem_results=temp("datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.scenicplus.dem.hdf5"),
        annotation_direct_path=temp("datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.scenicplus.annotation_direct.h5ad"),
        annotation_extended_path=temp("datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.scenicplus.annotation_extended.h5ad"),
        tf_names_path=temp("datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.scenicplus.tf_names.txt"),
        out="datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.scenicplus.tfb.csv"
    resources:
        mem_mb=512000,
        runtime=360,
    singularity:
        'workflow/envs/scenicplus.sif'
    shell:
        """
        python workflow/scripts/methods/scenicplus/tfb.py \
        -i {input.pre} \
        -p {input.p2g} \
        -d {input.mdata} \
        -r {input.cistarget_rankings_human} \
        -s {input.cistarget_scores_human} \
        -g {params.organism} \
        -t {output.cistarget_results} \
        -u {output.dem_results} \
        -o {output.out} \
        -c {params.n_cores} \
        --annotation_direct_path {output.annotation_direct_path}\
        --annotation_extended_path {output.annotation_extended_path}\
        --tf_names_path {output.tf_names_path}\
        --path_to_motif_annotations_human {input.path_to_motif_annotations_human}\
        --path_to_motif_annotations_mouse {input.path_to_motif_annotations_mouse}\
        --temp_dir {params.temp_dir} \
        --ray_tmp_dir {params.ray_tmp_dir}
        """


rule mdl_scenicplus:
    threads: 32
    input:
        pre=lambda wildcards: map_rules('pre', wildcards.pre),
        p2g=lambda wildcards: map_rules('p2g', wildcards.p2g),
        tfb=lambda wildcards: map_rules('tfb', wildcards.tfb),
        cistarget_rankings_human=rules.download_cistarget.output.human_rankings,
        cistarget_rankings_mouse=rules.download_cistarget.output.mouse_rankings
    output:
        out="datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.scenicplus.mdl.csv"
    params:
        tmp_dir="/tmp",
        n_cores=32,
        organism=lambda w: config['datasets'][w.dataset]['organism'],
        method_mdl=config['methods']['scenicplus']['method_mdl'],
        order_regions_to_genes_by=config['methods']['scenicplus']['order_regions_to_genes_by'],
        order_TFs_to_genes_by=config['methods']['scenicplus']['order_TFs_to_genes_by'],
        gsea_n_perm=config['methods']['scenicplus']['gsea_n_perm'],
        quantile_thresholds_region_to_gene=config['methods']['scenicplus']['quantile_thresholds_region_to_gene'],
        top_n_regionTogenes_per_gene=config['methods']['scenicplus']['top_n_regionTogenes_per_gene'],
        top_n_regionTogenes_per_region=config['methods']['scenicplus']['top_n_regionTogenes_per_region'],
        min_regions_per_gene=config['methods']['scenicplus']['min_regions_per_gene'],
        rho_threshold=config['methods']['scenicplus']['rho_threshold'],
        min_target_genes=config['methods']['scenicplus']['min_target_genes'],
    singularity:
        'workflow/envs/scenicplus.sif'
    shell:
        """
        python workflow/scripts/methods/scenicplus/mdl.py \
        -i {input.pre} \
        -p {input.p2g} \
        -t {input.tfb} \
        -l {input.cistarget_rankings_human} \
        -k {input.cistarget_rankings_mouse} \
        -o {output.out} \
        -m {params.method_mdl} \
        -c {params.n_cores} \
        -g {params.organism} \
        --temp_dir {params.tmp_dir} \
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


rule mdl_o_scenicplus:
    threads: 32
    input:
        # data
        frags=list_frags_files,
        mdata=rules.extract_case.output.mdata,
        # annotations
        chrom_sizes_h=rules.download_genomesizes.output.hg38,
        chrom_sizes_m=rules.download_genomesizes.output.mm10,
        annot_m=rules.download_gene_annotations.output.m,
        annot_h=rules.download_gene_annotations.output.h,
        path_to_motif_annotations_human=rules.download_motifs.output.h,
        path_to_motif_annotations_mouse=rules.download_motifs.output.m,
        # cistarget
        cistarget_rankings_human=rules.download_cistarget.output.human_rankings,
        cistarget_scores_human=rules.download_cistarget.output.human_scores,
    params:
        # general params
        n_cores=32,
        organism=lambda w: config['datasets'][w.dataset]['organism'],
        ray_tmp_dir="/tmp",
        temp_dir=temp(directory(local('datasets/{dataset}/cases/{case}/runs/scenicplus_tmp/src/'))),
        tmp_scenicplus=temp(directory(local('datasets/{dataset}/cases/{case}/runs/scenicplus_tmp'))), 
        # p2g params
        search_space_upstream="1000 150000",
        search_space_downstream="1000 150000",
        search_space_extend_tss="10 10",
        remove_promoters=False,  # Fixed to False in the snakemake pipeline
        use_gene_boundaries=False,  # Fixed to False in the snakemake pipeline
        region_to_gene_importance_method="GBM",
        region_to_gene_correlation_method="SR",
        # tfb params
        # mdl params
        method_mdl="GBM",
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
    output:
        out='datasets/{dataset}/cases/{case}/runs/o_scenicplus.o_scenicplus.o_scenicplus.o_scenicplus.mdl.csv',
        cistarget_results=temp("datasets/{dataset}/cases/{case}/runs/o_scenicplus.cistarget.hdf5"),
        dem_results=temp("datasets/{dataset}/cases/{case}/runs/o_scenicplus.dem.hdf5"),
        annotation_direct_path=temp("datasets/{dataset}/cases/{case}/runs/o_scenicplus.annotation_direct.h5ad"),
        annotation_extended_path=temp("datasets/{dataset}/cases/{case}/runs/o_scenicplus.annotation_extended.h5ad"),
        tf_names_path=temp("datasets/{dataset}/cases/{case}/runs/o_scenicplus.tf_names.txt"),
    shell:
        """
        python workflow/scripts/methods/scenicplus/src.py \
        -f {input.frags} \
        -i {input.mdata} \
        -m {input.chrom_sizes_m} \
        -j {input.chrom_sizes_h} \
        -o {output.out} \
        -t {params.temp_dir} \
        -g {params.organism} \
        -c {params.n_cores} \
        -a {input.annot_m} \
        -b {input.annot_h} \
        -r {input.cistarget_rankings_human} \
        -s {input.cistarget_scores_human} \
        --path_to_motif_annotations_human {input.path_to_motif_annotations_human}\
        --path_to_motif_annotations_mouse {input.path_to_motif_annotations_mouse}\
        --ray_tmp_dir {params.ray_tmp_dir} \
        --tmp_scenicplus {params.tmp_scenicplus} \
        --search_space_upstream "{params.search_space_upstream}" \
        --search_space_downstream "{params.search_space_downstream}" \
        --search_space_extend_tss "{params.search_space_extend_tss}" \
        --remove_promoters {params.remove_promoters} \
        --use_gene_boundaries {params.use_gene_boundaries} \
        --region_to_gene_importance_method {params.region_to_gene_importance_method} \
        --region_to_gene_correlation_method {params.region_to_gene_correlation_method} \
        --method_mdl {params.method_mdl} \
        --order_regions_to_genes_by {params.order_regions_to_genes_by} \
        --order_TFs_to_genes_by {params.order_TFs_to_genes_by} \
        --gsea_n_perm {params.gsea_n_perm} \
        --quantile_thresholds_region_to_gene {params.quantile_thresholds_region_to_gene} \
        --top_n_regionTogenes_per_gene {params.top_n_regionTogenes_per_gene} \
        --top_n_regionTogenes_per_region {params.top_n_regionTogenes_per_region} \
        --min_regions_per_gene {params.min_regions_per_gene} \
        --rho_threshold {params.rho_threshold} \
        --min_target_genes {params.min_target_genes} \
        --cistarget_results {output.cistarget_results} \
        --dem_results {output.dem_results} \
        --annotation_direct_path {output.annotation_direct_path}\
        --annotation_extended_path {output.annotation_extended_path}\
        --tf_names_path {output.tf_names_path}\
        """
        

