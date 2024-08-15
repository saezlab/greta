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
        k=config['methods']['celloracle']['k']
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
        thr_coaccess=config['methods']['celloracle']['thr_coaccess'],
        ext=config['methods']['celloracle']['ext']
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
        fpr=config['methods']['celloracle']['fpr'],
        blen=config['methods']['celloracle']['blen'],
        tfb_thr=config['methods']['celloracle']['tfb_thr']
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
        a=config['methods']['celloracle']['a'],
        p=config['methods']['celloracle']['p'],
        n=config['methods']['celloracle']['n'],
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

rule src_celloracle:
    input:
        'datasets/{dataset}/cases/{case}/mdata.h5mu',
    singularity:
        'workflow/envs/celloracle.sif'
    benchmark:
        'benchmarks/{dataset}.{case}.celloracle.src.txt',
    output:
        pp=temp(local('datasets/{dataset}/cases/{case}/runs/celloracle.src.peaks.csv')),
        pc=temp(local('datasets/{dataset}/cases/{case}/runs/celloracle.src.conns.csv')),
        gr='datasets/{dataset}/cases/{case}/runs/celloracle.src.csv',
    params:
        organism=lambda w: config['datasets'][w.dataset]['organism'],
        k=config['methods']['celloracle']['k']
        thr_coaccess=config['methods']['celloracle']['thr_coaccess'],
        ext=config['methods']['celloracle']['ext'],
        fpr=config['methods']['celloracle']['fpr'],
        blen=config['methods']['celloracle']['blen'],
        tfb_thr=config['methods']['celloracle']['tfb_thr'],
        a=config['methods']['celloracle']['a'],
        p=config['methods']['celloracle']['p'],
        n=config['methods']['celloracle']['n'],
    resources:
        mem_mb=256000,
        runtime=720,
    shell:
        """
        Rscript workflow/scripts/methods/celloracle/src.R \
        {input} \
        {params.organism} \
        {params.ext} \
        {output.pp} \
        {output.pc}

        python workflow/scripts/methods/celloracle/src.py \
        -a {input} \
        -b {output.pp} \
        -c {output.pc} \
        -d {params.organism} \
        -e {params.thr_coaccess} \
        -f {params.fpr} \
        -g {params.blen} \
        -i {params.tfb_thr} \
        -j {params.a} \
        -k {params.p} \
        -l {params.n} \
        -m {params.k} \
        -n {output.gr}
        """
