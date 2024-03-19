rule download_tfbs:
    output:
        h=directory('gdata/tfbs/hg38'),
        m=directory('gdata/tfbs/mm10'),
    shell:
        """
        wget 'https://www.embl.de/download/zaugg/diffTF/TFBS/TFBS_hg38_PWMScan_HOCOMOCOv11.tar.gz' -O {output.h}.tar.gz
        wget 'https://www.embl.de/download/zaugg/diffTF/TFBS/TFBS_mm10_PWMScan_HOCOMOCOv10.tar.gz' -O {output.m}.tar.gz
        mkdir {output.h} {output.m}
        tar -xvf {output.h}.tar.gz -C {output.h}
        tar -xvf {output.m}.tar.gz -C {output.m}
        rm {output.h}.tar.gz {output.m}.tar.gz
        mv {output.h}/*/* {output.h}
        mv {output.m}/*/* {output.m}
        rm -r {output.h}/PWMScan_* {output.m}/PWMScan_*
        """

rule pre_granie:
    input:
        'datasets/{dataset}/cases/{case}/mdata.h5mu'
    singularity:
        'workflow/envs/granie.sif'
    benchmark:
        'benchmarks/{dataset}.{case}.granie.pre.txt'
    output:
        'datasets/{dataset}/cases/{case}/runs/granie.pre.h5mu'
    shell:
        """
        python workflow/scripts/methods/granie/pre.py \
        -i {input} \
        -o {output}
        Rscript workflow/scripts/methods/granie/pre.R \
        {output}
        """

rule p2g_granie:
    input:
        d='datasets/{dataset}/cases/{case}/runs/{pre}.pre.h5mu',
        g='gdata/geneids',
    singularity:
        'workflow/envs/granie.sif'
    benchmark:
        'benchmarks/{dataset}.{case}.{pre}.granie.p2g.txt'
    output:
        t=temp(directory(local('datasets/{dataset}/cases/{case}/runs/{pre}.granie_tmp'))),
        p='datasets/{dataset}/cases/{case}/runs/{pre}.granie.p2g.csv'
    params:
        organism=lambda w: config['datasets'][w.dataset]['organism'],
        ext=500000,
        thr_fdr=0.2, # Need to change back to 0.2
    shell:
        """
        Rscript workflow/scripts/methods/granie/p2g.R \
        {input.d} \
        {params.organism} \
        {input.g} \
        {output.t} \
        {params.ext} \
        {output.p}
        """

rule tfb_granie:
    input:
        d='datasets/{dataset}/cases/{case}/runs/{pre}.pre.h5mu',
        p='datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.p2g.csv',
        g='gdata/geneids',
        t='gdata/tfbs',
    singularity:
        'workflow/envs/granie.sif'
    benchmark:
        'benchmarks/{dataset}.{case}.{pre}.{p2g}.granie.tfb.txt'
    output:
        t=temp(directory(local('datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.granie_tmp'))),
        p='datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.granie.tfb.csv'
    params:
        organism=lambda w: config['datasets'][w.dataset]['organism'],
    shell:
        """
        Rscript workflow/scripts/methods/granie/tfb.R \
        {input.d} \
        {params.organism} \
        {input.g} \
        {output.t} \
        {input.t} \
        {input.p} \
        {output.p}
        """