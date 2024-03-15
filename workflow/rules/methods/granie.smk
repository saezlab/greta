rule download_tfbs:
    output:
        h=directory('gdata/tfbs/TFBS_hg38_PWMScan_HOCOMOCOv11'),
        m=directory('gdata/tfbs/TFBS_mm10_PWMScan_HOCOMOCOv10'),
    shell:
        """
        wget 'https://www.embl.de/download/zaugg/diffTF/TFBS/TFBS_hg38_PWMScan_HOCOMOCOv11.tar.gz' -O {output.h}.tar.gz
        wget 'https://www.embl.de/download/zaugg/diffTF/TFBS/TFBS_mm10_PWMScan_HOCOMOCOv10.tar.gz' -O {output.m}.tar.gz
        mkdir {output.h} {output.m}
        tar -xvf {output.h}.tar.gz -C {output.h}
        tar -xvf {output.m}.tar.gz -C {output.m}
        rm {output.h}.tar.gz {output.m}.tar.gz
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
        thr_fdr=0.5, # Need to change back to 0.2
    shell:
        """
        Rscript workflow/scripts/methods/granie/p2g.R \
        {input.d} \
        {params.organism} \
        {input.g} \
        {output.t} \
        {params.ext} \
        {params.thr_fdr} \
        {output.p}
        """
