rule extract_tss_figr:
    output:
        "gdata/alltss/figr.csv"
    singulatiry:
	'workflow/envs/figr.sif'
    shell:
        """
        Rscript scripts/analysis/tss/figr.R {output}
        """

rule extract_tss_dictys:
    conda:
        '../../envs/dictys.yaml'
    params:
        url_gtf="http://ftp.ensembl.org/pub/release-107/gtf/homo_sapiens/Homo_sapiens.GRCh38.107.gtf.gz"
    output:
        gtf=temp(local('gdata/alltss/gene.gtf')),
        i=temp(local('gdata/alltss/gene.bed')),
        o='gdata/alltss/dictys.csv'
    shell:
        """
        wget -O {output.gtf}.gz {params.url_gtf} && \
        gunzip {output.gtf}.gz
        dictys_helper gene_gtf.sh {output.gtf} {output.i}
        # Further preprocessing
        python3 scripts/analysis/tss/dictys.py \
        -i {output.i} \
        -o {output.o}
        """

rule extract_tss_celloracle:
    output:
        o="gdata/alltss/celloracle.csv"
    singulatiry:
        'workflow/envs/celloracle.sif'
    shell:
        """
        python scripts/analysis/tss/celloracle.py \
        -o {output.o}
        """

rule extract_tss_hummus:
    output:
        "gdata/alltss/hummus.csv"
    singulatiry:
        'workflow/envs/hummus.sif'
    shell:
        """
        Rscript scripts/analysis/tss/hummus.R {output}
        """

rule extract_tss_pando:
    output:
        "gdata/alltss/pando.csv"
    singulatiry:
        'workflow/envs/pando.sif'
    shell:
        """
        Rscript scripts/analysis/tss/pando.R {output}
        """

rule extract_tss_granie:
    output:
        "gdata/alltss/granie.csv"
    singulatiry:
        'workflow/envs/granie.sif'
    shell:
        """
        Rscript scripts/analysis/tss/granie.R {output}
        """

rule extract_tss_scenicplus:
    output:
        h="gdata/alltss/scenicplus.csv"
    singulatiry:
        'workflow/envs/scenicplus.sif'
    shell:
        """
        python scripts/analysis/tss/scenicplus.py \
        -j {output.h}
        """

rule compare_tss:
    input:
        g="gdata/alltss/granie.csv",
        s="gdata/alltss/scenicplus.csv",
        c="gdata/alltss/celloracle.csv",
        m="gdata/alltss/hummus.csv",
	p="gdata/alltss/pando.csv",
        f="gdata/alltss/figr.csv",
	d="gdata/alltss/dictys.csv"
    output:
        o="/analysis/alltss/ocoef.csv"
    shell:
        """
        # Process datasets (filtering). Getting Overlap Coefficient
        python scripts/analysis/tss/tss_result.py \
        -g {input.g} \
        -s {input.s} \
        -c {input.c} \
        -p {input.p} \
        -f {input.f} \
	-m {input.m} \
	-d {input.d} \
        -o {output.o} 
        echo "Done"
        """

