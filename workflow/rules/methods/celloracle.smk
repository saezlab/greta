rule pre_celloracle:
    input:
        'datasets/{dataset}/cases/{case}/mdata.h5mu'
    singularity:
        'workflow/envs/celloracle.sif'
    benchmark:
        'benchmarks/{dataset}.{case}.celloracle.{p2g}.{tfb}.{mdl}.pre.txt'
    output:
        'datasets/{dataset}/cases/{case}/runs/celloracle.{p2g}.{tfb}.{mdl}.pre.h5mu'
    params:
        k=20
    shell:
        """
        python workflow/scripts/methods/celloracle/pre.py \
        -i {input} \
        -k {params.k} \
        -o {output}
        """

rule download_genomesizes:
    output:
        hg38='gdata/sizes/hg38.txt',
        mm10='gdata/sizes/mm10.txt',
        size=directory('gdata/sizes')
    shell:
        """
        wget 'http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes' -O {output.hg38}
        wget 'http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes' -O {output.mm10}
        """

rule p2g_celloracle:
    input:
        data='datasets/{dataset}/cases/{case}/runs/{pre}.celloracle.{tfb}.{mdl}.pre.h5mu',
        sizes='gdata/sizes/',
    singularity:
        'workflow/envs/celloracle.sif'
    benchmark:
        'benchmarks/{dataset}.{case}.{pre}.celloracle.{tfb}.{mdl}.p2g.txt'
    output:
        pp='datasets/{dataset}/cases/{case}/runs/{pre}.celloracle.{tfb}.{mdl}.all_peaks.csv',
        pc='datasets/{dataset}/cases/{case}/runs/{pre}.celloracle.{tfb}.{mdl}.cicero_connections.csv',
        pg='datasets/{dataset}/cases/{case}/runs/{pre}.celloracle.{tfb}.{mdl}.p2g.csv',
    params:
        organism=lambda w: config['datasets'][w.dataset]['organism'],
        k=50,
        thr_coaccess=0.8,
    shell:
        """
        Rscript workflow/scripts/methods/celloracle/p2g.R \
        {input.data} \
        {params.organism} \
        {params.k} \
        {output.pp} \
        {output.pc}

        python workflow/scripts/methods/celloracle/p2g.py \
        -a {output.pp} \
        -c {output.pc} \
        -o {params.organism} \
        -t {params.thr_coaccess} \
        -p {output.pg}
        """

rule download_genomes:
    output:
        d=directory('gdata/genomes'),
    singularity:
        'workflow/envs/celloracle.sif'
    shell:
        """
        python workflow/scripts/methods/celloracle/download_genomes.py \
        -d {output.d}
        """

rule tfb_celloracle:
    input:
        p2g='datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.celloracle.{mdl}.p2g.csv',
        gdir='gdata/genomes',
    singularity:
        'workflow/envs/celloracle.sif'
    benchmark:
        'benchmarks/{dataset}.{case}.{pre}.{p2g}.celloracle.{mdl}.tfb.txt',
    output:
        'datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.celloracle.{mdl}.tfb.csv'
    params:
        organism=lambda w: config['datasets'][w.dataset]['organism'],
        fpr=0.02,
        blen=200,
        tfb_thr=10,
    shell:
        """
        python workflow/scripts/methods/celloracle/tfb.py \
        -p {input.p2g} \
        -g {params.organism} \
        -f {params.fpr} \
        -b {params.blen} \
        -t {params.tfb_thr} \
        -o {output}
        """

rule mdl_celloracle:
    input:
        pre='datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.celloracle.pre.h5mu',
        p2g='datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.celloracle.p2g.csv',
        tfb='datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.celloracle.tfb.csv',
    singularity:
        'workflow/envs/celloracle.sif'
    benchmark:
        'benchmarks/{dataset}.{case}.{pre}.{p2g}.{tfb}.celloracle.mdl.txt',
    output:
        'datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.celloracle.grn.csv'
    params:
        a=10,
        p=0.001,
        n=2000,
    shell:
        """
        python workflow/scripts/methods/celloracle/mdl.py \
        -m {input.pre} \
        -g {input.p2g} \
        -t {input.tfb} \
        -a {params.a} \
        -p {params.p} \
        -n {params.n} \
        -o {output}
        """
