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
        

rule pre_scenicplus:
    input:
        frags = 'datasets/{dataset}/smpl.frags.tsv.gz',
        mudata = 'datasets/{dataset}/cases/{case}/mdata.h5mu'
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
        -o {output.d} \
        -t {output.tmp_scenicplus} \
        -g {params.organism} \
        -n {params.n_cores}
        """


# motif_enrichment_cistarget
# download_genome_annotations
# motif_enrichment_dem, prepare_menr, get_search_space, region_to_gene  
# tf_to_gene, eGRN_direct, eGRN_extended
