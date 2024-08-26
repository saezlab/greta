rule download_motifs:
    input:
        url_h = "https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl",
        url_m = "https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl"
    output:
        h = "aertslab/motifs-v10nr_clust/nr.hgnc-m0.001-o0.0.tbl",
        m = "aertslab/motifs-v10nr_clust/nr.mgi-m0.001-o0.0.tbl"
    shell:
        """
        wget -O {output.h} {input.url_h}
        wget -O {output.m} {input.url_m}
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
        human_rankings = "aertslab/cistarget/human_motif_SCREEN_rankings.feather",
        human_scores = "aertslab/cistarget/human_motif_SCREEN_scores.feather",
        mouse_rankings = "aertslab/cistarget/mouse_motif_SCREEN_rankings.feather",
        mouse_scores = "aertslab/cistarget/mouse_motif_SCREEN_scores.feather"
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
        tmp_scenicplus = temp(directory(local('datasets/{dataset}/cases/{case}/runs/scenicplus_tmp')))
    singularity:
        'workflow/envs/scenicplus.sif'
    params:
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
        -t {output.tmp_scenicplus} \
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
        temp_dir = temp(directory(local('datasets/{dataset}/cases/{case}/runs/{pre}.scenicplus_tmp'))),
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
        pre = "datasets/{dataset}/cases/{case}/runs/{pre}.pre.h5mu",
        p2g = "datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.p2g.csv",
        cistarget_rankings_human = "aertslab/cistarget/human_motif_SCREEN_rankings.feather",
        cistarget_scores_human = "aertslab/cistarget/human_motif_SCREEN_scores.feather",
    params:
        organism=lambda w: config['datasets'][w.dataset]['organism'],
    output:
        tfb = "datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.scenicplus.tfb.csv"
    singularity:
        'workflow/envs/scenicplus.sif'
    shell:
        """
        python workflow/scripts/methods/scenicplus/tfb.py \
        -i {input.pre} \
        -p {input.p2g} \
        -r {input.cistarget_rankings_human} \
        -s {input.cistarget_scores_human} \
        -g {params.organism} \
        -o {output.tfb}
        """
# motif_enrichment_cistarget
# download_genome_annotations
# motif_enrichment_dem, prepare_menr, get_search_space, region_to_gene  
# tf_to_gene, eGRN_direct, eGRN_extended
