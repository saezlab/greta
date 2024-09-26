localrules: download_tfbs


rule download_tfbs:
    output:
        h=directory('gdata/tfbs/hg38'),
        m=directory('gdata/tfbs/mm10'),
        d=directory('gdata/tfbs'),
    shell:
        """
        wget 'https://s3.embl.de/zaugg-web/GRaNIE/TFBS/hg38/PWMScan_HOCOMOCOv12_H12INVIVO.tar.gz' -O {output.h}.tar.gz
        wget 'https://s3.embl.de/zaugg-web/GRaNIE/TFBS/mm10/PWMScan_HOCOMOCOv12_H12INVIVO.tar.gz' -O {output.m}.tar.gz
        mkdir {output.h} {output.m}
        tar -xvf {output.h}.tar.gz -C {output.h}
        tar -xvf {output.m}.tar.gz -C {output.m}
        rm {output.h}.tar.gz {output.m}.tar.gz
        mv {output.h}/*/*/* {output.h}
        mv {output.m}/*/* {output.m}
        rm -r {output.h}/PWMScan_* {output.m}/pwmscan_*
        """


rule pre_granie:
    threads: 32
    input:
        mdata=rules.extract_case.output.mdata,
    singularity:
        'workflow/envs/granie.sif'
    output:
        tmp=temp(local('datasets/{dataset}/cases/{case}/runs/granie.pre_tmp.h5mu')),
        out='datasets/{dataset}/cases/{case}/runs/granie.pre.h5mu'
    shell:
        """
        python workflow/scripts/methods/granie/pre.py \
        -i {input.mdata} \
        -o {output.tmp}
        Rscript workflow/scripts/methods/granie/pre.R \
        {output.tmp}
        python workflow/scripts/methods/granie/pre_post.py \
        -i {output.tmp} \
        -o {output.out}
        """


rule p2g_granie:
    threads: 32
    input:
        pre=lambda wildcards: map_rules('pre', wildcards.pre),
        g=rules.download_geneids.output.dr,
    singularity:
        'workflow/envs/granie.sif'
    output:
        t=temp(directory(local('datasets/{dataset}/cases/{case}/runs/{pre}.granie_tmp'))),
        out='datasets/{dataset}/cases/{case}/runs/{pre}.granie.p2g.csv'
    params:
        organism=lambda w: config['datasets'][w.dataset]['organism'],
        ext=config['methods']['granie']['ext'],
    resources:
        mem_mb=128000,
        runtime=1440,
    shell:
        """
        Rscript workflow/scripts/methods/granie/p2g.R \
        {input.pre} \
        {params.organism} \
        {input.g} \
        {output.t} \
        {params.ext} \
        {output.out}
        """


rule tfb_granie:
    threads: 32
    input:
        pre=lambda wildcards: map_rules('pre', wildcards.pre),
        p2g=lambda wildcards: map_rules('p2g', wildcards.p2g),
        g=rules.download_geneids.output.dr,
        t=rules.download_tfbs.output.d,
    singularity:
        'workflow/envs/granie.sif'
    output:
        t=temp(directory(local('datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.granie_tmp'))),
        out='datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.granie.tfb.csv'
    params:
        organism=lambda w: config['datasets'][w.dataset]['organism'],
    resources:
        mem_mb=128000,
    shell:
        """
        Rscript workflow/scripts/methods/granie/tfb.R \
        {input.pre} \
        {params.organism} \
        {input.g} \
        {output.t} \
        {input.t} \
        {input.p2g} \
        {output.out}
        """


rule mdl_granie:
    threads: 32
    input:
        pre=lambda wildcards: map_rules('pre', wildcards.pre),
        p2g=lambda wildcards: map_rules('p2g', wildcards.p2g),
        tfb=lambda wildcards: map_rules('tfb', wildcards.tfb),
        g=rules.download_geneids.output.dr,
    singularity:
        'workflow/envs/granie.sif'
    output:
        t=temp(directory(local('datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.granie_tmp'))),
        out='datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.granie.mdl.csv'
    params:
        organism=lambda w: config['datasets'][w.dataset]['organism'],
        thr_fdr=config['methods']['granie']['thr_fdr'],
    resources:
        mem_mb=128000,
        runtime=720,
    shell:
        """
        Rscript workflow/scripts/methods/granie/mdl.R \
        {input.pre} \
        {params.organism} \
        {input.g} \
        {output.t} \
        {input.p2g} \
        {input.tfb} \
        {params.thr_fdr} \
        {output.out}
        """


rule mdl_o_granie:
    threads: 32
    input:
        mdata=rules.extract_case.output.mdata,
        g=rules.download_geneids.output.dr,
        t=rules.download_tfbs.output.d,
    singularity:
        'workflow/envs/granie.sif'
    output:
        h=temp(local('datasets/{dataset}/cases/{case}/runs/pre.granie.src.h5mu')),
        t=temp(directory(local('datasets/{dataset}/cases/{case}/runs/granie_tmp.src'))),
        out='datasets/{dataset}/cases/{case}/runs/o_granie.o_granie.o_granie.o_granie.grn.csv'
    params:
        organism=lambda w: config['datasets'][w.dataset]['organism'],
        ext=config['methods']['granie']['ext'],
        thr_fdr=config['methods']['granie']['thr_fdr'],
    resources:
        mem_mb=128000,
    shell:
        """
        python workflow/scripts/methods/granie/pre.py \
        -i {input.mdata} \
        -o {output.h}

        Rscript workflow/scripts/methods/granie/src.R \
        {output.h} \
        {params.organism} \
        {input.g} \
        {output.t} \
        {params.ext} \
        {input.t} \
        {params.thr_fdr} \
        {output.out}
        """
