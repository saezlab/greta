localrules: download_motifs, download_gene_annotations_scenicplus


rule download_motifs:
    threads: 1
    params:
        url_h="https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl",
        url_m="https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl"
    output:
        h="gdata/aertslab/motifs-v10nr_clust/nr.hgnc-m0.001-o0.0.tbl",
        m="gdata/aertslab/motifs-v10nr_clust/nr.mgi-m0.001-o0.0.tbl"
    shell:
        """
        wget -O {output.h} {params.url_h}
        wget -O {output.m} {params.url_m}
        """


rule download_gene_annotations_scenicplus:
    threads: 1
    output:
        gann=expand('gdata/aertslab/genomes/{org}/{org}_gannot.bed', org=['hg38', 'mm10']),
        cist=expand('gdata/aertslab/genomes/{org}/cist_{org}_gannot.bed', org=['hg38', 'mm10']),
        csiz=expand('gdata/aertslab/genomes/{org}/{org}_csizes.tsv', org=['hg38', 'mm10']),
    singularity:
        "workflow/envs/scenicplus.sif"
    shell:
        """
        python workflow/scripts/methods/scenicplus/download_gene_annot.py \
        -g {output.gann} \
        -c {output.cist} \
        -s {output.csiz}
        """


rule download_cistarget:
    threads: 1
    output:
        human_rankings="gdata/aertslab/cistarget/human_motif_SCREEN.regions_vs_motifs.rankings.feather",
        human_scores="gdata/aertslab/cistarget/human_motif_SCREEN.regions_vs_motifs.scores.feather",
        mouse_rankings="gdata/aertslab/cistarget/mouse_motif_SCREEN.regions_vs_motifs.rankings.feather",
        mouse_scores="gdata/aertslab/cistarget/mouse_motif_SCREEN.regions_vs_motifs.scores.feather"
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
        cist=rules.download_gene_annotations_scenicplus.output.cist,
        csiz=rules.download_gene_annotations_scenicplus.output.csiz,
    output:
        tmp_scenicplus=directory(local('datasets/{dataset}/cases/{case}/runs/scenicplus_tmp')),
        out='datasets/{dataset}/cases/{case}/runs/scenicplus.pre.h5mu',
        cistopic_obj='datasets/{dataset}/cases/{case}/runs/citopic.pkl'
    singularity:
        'workflow/envs/scenicplus.sif'
    params:
        ray_tmp_dir="/tmp",
        organism=lambda w: config['datasets'][w.dataset]['organism'],
    resources:
        mem_mb=256000,
        runtime=config['max_mins_per_step'],
    shell:
        """
        python workflow/scripts/methods/scenicplus/pre.py \
        -f {input.frags} \
        -i {input.mdata} \
        -o {output.out} \
        -t {output.tmp_scenicplus} \
        -g {params.organism} \
        -n {threads} \
        -a {input.cist} \
        -b {input.csiz} \
        --ray_tmp_dir '{params.ray_tmp_dir}' \
        -c {output.cistopic_obj}
        """


rule p2g_scenicplus:
    threads: 32
    input:
        pre=lambda wildcards: map_rules('pre', wildcards.pre),
        csz=rules.gen_genome_celloracle.output,
    singularity:
        'workflow/envs/scenicplus.sif'
    params:
        temp_dir=temp(directory(local('datasets/{dataset}/cases/{case}/runs/scenicplus_tmp/{pre}'))),
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
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    shell:
        """
        python workflow/scripts/methods/scenicplus/p2g.py \
        -i {input.pre} \
        -c {threads} \
        -t {params.temp_dir} \
        -m {input.csz} \
        -j {input.csz} \
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
        annot_m=rules.download_gene_annotations_scenicplus.output.gann,
        annot_h=rules.download_gene_annotations_scenicplus.output.gann,
    params:
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
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
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
        -c {threads} \
        --annotation_direct_path {output.annotation_direct_path}\
        --annotation_extended_path {output.annotation_extended_path}\
        --tf_names_path {output.tf_names_path}\
        --path_to_motif_annotations_human {input.path_to_motif_annotations_human}\
        --path_to_motif_annotations_mouse {input.path_to_motif_annotations_mouse}\
        --temp_dir {params.temp_dir} \
        --ray_tmp_dir {params.ray_tmp_dir} \
        --gannot {input.annot_h}
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
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
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
        -c {threads} \
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
        mdata='datasets/{dataset}/cases/{case}/mdata.h5mu',
#        mdata=rules.extract_case.output.mdata,
        # annotations
        csz=rules.gen_genome_celloracle.output,
#        annot_m=rules.download_gene_annotations.output.m,
#        annot_h=rules.download_gene_annotations.output.h,
        path_to_motif_annotations_human=rules.download_motifs.output.h,
        path_to_motif_annotations_mouse=rules.download_motifs.output.m,
        # cistarget
        cistarget_rankings_human=rules.download_cistarget.output.human_rankings,
        cistarget_scores_human=rules.download_cistarget.output.human_scores,
        # cistopic
        cistopic_obj=rules.pre_scenicplus.output.cistopic_obj
    params:
        organism=lambda w: config['datasets'][w.dataset]['organism'],
        temp_numba_cache_dir=".scenicplus_tmp",
        ray_tmp_dir="/tmp",
        tmp_dir=directory(local('datasets/{dataset}/cases/{case}/runs/scenicplus_tmp')),
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
        # Only download in scenicplus snakemake if not already existing
        annot_m = "gdata/aertslab/genomes/mm10/o_scenicplus_mm10_ensdb_v79.tsv",
        annot_h = "gdata/aertslab/genomes/hg38/o_scenicplus_hg38_ensdb_v86.tsv",
        chrom_sizes_h="gdata/aertslab/genomes/hg38/o_scenicplus_hg38_chromsizes.tsv",
        chrom_sizes_m="gdata/aertslab/genomes/mm10/o_scenicplus_mm10_chromsizes.tsv",
        cistarget_results=temp("datasets/{dataset}/cases/{case}/runs/o_scenicplus.cistarget.hdf5"),
        dem_results=temp("datasets/{dataset}/cases/{case}/runs/o_scenicplus.dem.hdf5"),
        annotation_direct_path=temp("datasets/{dataset}/cases/{case}/runs/o_scenicplus.annotation_direct.h5ad"),
        annotation_extended_path=temp("datasets/{dataset}/cases/{case}/runs/o_scenicplus.annotation_extended.h5ad"),
        tf_names_path=temp("datasets/{dataset}/cases/{case}/runs/o_scenicplus.tf_names.txt"),
        search_space_path = temp("datasets/{dataset}/cases/{case}/runs/o_scenicplus.search_space.tsv"),
    singularity:
        'workflow/envs/scenicplus.sif'
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    output:
        out='datasets/{dataset}/cases/{case}/runs/o_scenicplus.o_scenicplus.o_scenicplus.o_scenicplus.mdl.csv',
#        cistarget_results=temp("datasets/{dataset}/cases/{case}/runs/o_scenicplus.cistarget.hdf5"),
#        dem_results=temp("datasets/{dataset}/cases/{case}/runs/o_scenicplus.dem.hdf5"),
#        annotation_direct_path=temp("datasets/{dataset}/cases/{case}/runs/o_scenicplus.annotation_direct.h5ad"),
#        annotation_extended_path=temp("datasets/{dataset}/cases/{case}/runs/o_scenicplus.annotation_extended.h5ad"),
#        tf_names_path=temp("datasets/{dataset}/cases/{case}/runs/o_scenicplus.tf_names.txt"),
#        search_space_path = temp("datasets/{dataset}/cases/{case}/runs/o_scenicplus.search_space.tsv"),
    shell:
        """
        export NUMBA_CACHE_DIR=$(pwd)/{params.tmp_dir}      

        scplus_pipeline=scplus_pipeline_{wildcards.dataset}_{wildcards.case}
        #mkdir -p  $scplus_pipeline
        rm -fr $scplus_pipeline
        scenicplus init_snakemake --out_dir $scplus_pipeline

        python workflow/scripts/methods/scenicplus/src_config.py \
        -f {input.frags} \
        -i {input.mdata} \
        -m {params.csz} \
        -j {params.csz} \
        --cistopic_path {input.cistopic_obj} \
        -o {output.out} \
        -t {params.tmp_dir} \
        -g {params.organism} \
        -c {threads} \
        -a {params.annot_m} \
        -b {params.annot_h} \
        -r {input.cistarget_rankings_human} \
        -s {input.cistarget_scores_human} \
        --path_to_motif_annotations_human {input.path_to_motif_annotations_human}\
        --path_to_motif_annotations_mouse {input.path_to_motif_annotations_mouse}\
        --ray_tmp_dir {params.ray_tmp_dir} \
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
        --cistarget_results {params.cistarget_results} \
        --dem_results {params.dem_results} \
        --annotation_direct_path {params.annotation_direct_path}\
        --annotation_extended_path {params.annotation_extended_path}\
        --tf_names_path {params.tf_names_path}\
        --search_space_path {params.search_space_path}\
        --output_config $scplus_pipeline/Snakemake/config/config.yaml

        cd $scplus_pipeline/Snakemake/
        snakemake --cores {threads} --keep-incomplete
        
        cd ../../
        python workflow/scripts/methods/scenicplus/src_merging.py \
        --grn_extended {params.tmp_dir}/o_scenicplus_eRegulons_extended.tsv\
        --grn_direct {params.tmp_dir}/o_scenicplus_eRegulons_direct.tsv\
        --output {output.out}
        
        rm -r $scplus_pipeline/
        """
