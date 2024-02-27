rule download_genomesizes:
    output:
        hg38='genomes/sizes/hg38.txt',
        mm10='genomes/sizes/mm10.txt',
        size=directory('genomes/sizes')
    shell:
        """
        wget 'http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes' -O {output.hg38}
        wget 'http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes' -O {output.mm10}
        """

rule celloracle_p2g_0:
    input:
        data='datasets/{dataset}/cases/{case}/mdata.h5mu',
        sizes='genomes/sizes/',
    singularity:
        'workflow/envs/celloracle.sif'
    benchmark:
        'benchmarks/celloracle/{dataset}.{case}.p2g.0.txt'
    output:
        path_all_peaks='datasets/{dataset}/cases/{case}/celloracle/all_peaks.csv',
        path_connections='datasets/{dataset}/cases/{case}/celloracle/cicero_connections.csv'
    params:
        organism=lambda w: config['datasets'][w.dataset]['organism'],
        k=50
    shell:
        """
        Rscript workflow/scripts/methods/celloracle/peak_corr.R \
        {input.data} \
        {params.organism} \
        {params.k} \
        {output.path_all_peaks} \
        {output.path_connections}
        """

rule celloracle_p2g_1:
    input:
        all_peaks='datasets/{dataset}/cases/{case}/celloracle/all_peaks.csv',
        connections='datasets/{dataset}/cases/{case}/celloracle/cicero_connections.csv'
    singularity:
        'workflow/envs/celloracle.sif'
    benchmark:
        'benchmarks/celloracle/{dataset}.{case}.p2g.1.txt'
    output:
        'datasets/{dataset}/cases/{case}/celloracle/p2g.csv'
    params:
        organism=lambda w: config['datasets'][w.dataset]['organism'],
        thr_coaccess=0.8
    shell:
        """
        python workflow/scripts/methods/celloracle/tss_annotation.py \
        -a {input.all_peaks} \
        -c {input.connections} \
        -o {params.organism} \
        -t {params.thr_coaccess} \
        -p {output}
        """

