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

