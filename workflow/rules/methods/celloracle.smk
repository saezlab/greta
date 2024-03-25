rule pre_celloracle:
    input:
        'datasets/{dataset}/cases/{case}/mdata.h5mu'
    singularity:
        'workflow/envs/celloracle.sif'
    benchmark:
        'benchmarks/{dataset}.{case}.celloracle.pre.txt'
    output:
        'datasets/{dataset}/cases/{case}/runs/celloracle.pre.h5mu'
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
        data='datasets/{dataset}/cases/{case}/runs/{pre}.pre.h5mu',
        sizes='gdata/sizes/',
    singularity:
        'workflow/envs/celloracle.sif'
    benchmark:
        'benchmarks/{dataset}.{case}.{pre}.celloracle.p2g.txt'
    output:
        pp=temp(local('datasets/{dataset}/cases/{case}/runs/{pre}.celloracle.peaks.csv')),
        pc=temp(local('datasets/{dataset}/cases/{case}/runs/{pre}.celloracle.conns.csv')),
        pg='datasets/{dataset}/cases/{case}/runs/{pre}.celloracle.p2g.csv',
    params:
        organism=lambda w: config['datasets'][w.dataset]['organism'],
        thr_coaccess=0.8,
        ext=500000
    shell:
        """
        Rscript workflow/scripts/methods/celloracle/p2g.R \
        {input.data} \
        {params.organism} \
        {params.ext} \
        {output.pp} \
        {output.pc}

        python workflow/scripts/methods/celloracle/p2g.py \
        -d {input.data} \
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
        d='datasets/{dataset}/cases/{case}/runs/{pre}.pre.h5mu',
        p='datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.p2g.csv',
        g='gdata/genomes',
    singularity:
        'workflow/envs/celloracle.sif'
    benchmark:
        'benchmarks/{dataset}.{case}.{pre}.{p2g}.celloracle.tfb.txt',
    output:
        'datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.celloracle.tfb.csv'
    params:
        organism=lambda w: config['datasets'][w.dataset]['organism'],
        fpr=0.02,
        blen=200,
        tfb_thr=10,
    shell:
        """
        python workflow/scripts/methods/celloracle/tfb.py \
        -d {input.d} \
        -p {input.p} \
        -g {params.organism} \
        -f {params.fpr} \
        -b {params.blen} \
        -t {params.tfb_thr} \
        -o {output}
        """

rule mdl_celloracle:
    input:
        pre='datasets/{dataset}/cases/{case}/runs/{pre}.pre.h5mu',
        p2g='datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.p2g.csv',
        tfb='datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.tfb.csv',
    singularity:
        'workflow/envs/celloracle.sif'
    benchmark:
        'benchmarks/{dataset}.{case}.{pre}.{p2g}.{tfb}.celloracle.mdl.txt',
    output:
        'datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.celloracle.mdl.csv'
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
