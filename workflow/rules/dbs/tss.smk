localrules: gen_tss_celloracle, gen_tss_crema, gen_tss_dictys, gen_tss_figr, gen_tss_granie, gen_tss_pando, gen_tss_scenicplus, gen_tss_scdori, gen_tss_scmtni, gen_tss_promoters


rule gen_tss_celloracle:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: 'workflow/envs/gretabench.sif'
    output: 'dbs/hg38/gen/tss/celloracle.bed.gz'
    shell:
        """
        python workflow/scripts/dbs/gen/tss/celloracle.py -o {output}
        """

rule gen_tss_crema:
    threads: 1
    singularity: 'workflow/envs/crema.sif'
    input: 'workflow/envs/crema.sif'
    output: 'dbs/hg38/gen/tss/crema.bed.gz'
    shell:
        """
        Rscript workflow/scripts/dbs/gen/tss/crema.R {output}
        """

rule gen_tss_dictys:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: rules.gen_ann_dictys.output
    output: 'dbs/hg38/gen/tss/dictys.bed.gz'
    shell:
        """
        python workflow/scripts/dbs/gen/tss/dictys.py \
        -i {input} \
        -o {output}
        """

rule gen_tss_directnet:
    threads: 1
    singularity: 'workflow/envs/directnet.sif'
    input: rules.gen_genome_inferelator.output.gtf
    output: 'dbs/hg38/gen/tss/directnet.bed.gz'
    shell:
        """
        python workflow/scripts/dbs/gen/tss/directnet.py {input} {output}
        """

rule gen_tss_figr:
    threads: 1
    singularity: 'workflow/envs/figr.sif'
    input: 'workflow/envs/figr.sif'
    output: 'dbs/hg38/gen/tss/figr.bed.gz'
    shell:
        """
        Rscript workflow/scripts/dbs/gen/tss/figr.R {output}
        """

rule gen_tss_granie:
    threads: 1
    singularity: 'workflow/envs/granie.sif'
    input: 'workflow/envs/granie.sif'
    output: 'dbs/hg38/gen/tss/granie.bed.gz'
    shell:
        """
        Rscript workflow/scripts/dbs/gen/tss/granie.R {output}
        """

rule gen_tss_hummus:
    threads: 1
    singularity: 'workflow/envs/hummus.sif'
    input: 'workflow/envs/hummus.sif'
    output: 'dbs/hg38/gen/tss/hummus.bed.gz'
    shell:
        """
        Rscript workflow/scripts/dbs/gen/tss/hummus.R {output}
        """

rule gen_tss_inferelator:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: rules.gen_tss_directnet.output
    output: 'dbs/hg38/gen/tss/inferelator.bed.gz'
    shell:
        """
        cp {input} {output}
        """

rule gen_tss_pando:
    threads: 1
    singularity: 'workflow/envs/pando.sif'
    input: 'workflow/envs/pando.sif'
    output: 'dbs/hg38/gen/tss/pando.bed.gz'
    shell:
        """
        Rscript workflow/scripts/dbs/gen/tss/pando.R {output}
        """

rule gen_tss_scdori:
    threads: 1
    conda: '../../envs/scdori.yaml'
    output: 'dbs/hg38/gen/tss/scdori.bed.gz'
    shell:
        """
        wget --no-verbose \
        'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.primary_assembly.basic.annotation.gtf.gz' \
        -O {output}.tmp
        python workflow/scripts/dbs/gen/tss/scdori.py {output}.tmp {output}
        rm {output}.tmp
        """

rule gen_tss_scenicplus:
    threads: 1
    singularity: 'workflow/envs/scenicplus.sif'
    input:
        img='workflow/envs/scenicplus.sif',
        tss=rules.gen_genome_scenicplus.output.tss,
    output: 'dbs/hg38/gen/tss/scenicplus.bed.gz'
    shell:
        """
        awk -v OFS='\\t' 'NR > 1 && $1 ~ /^chr/ && $4 != "." {{print $1, $2, $3, $4}}' {input.tss} | \
        gzip -c > {output}
        """

rule gen_tss_scmtni:
    threads: 1
    singularity: 'workflow/envs/scmtni.sif'
    input: rules.gen_motif_scmtni.output.promoters_dir,
    output: 'dbs/hg38/gen/tss/scmtni.bed.gz'
    shell:
        """
        python workflow/scripts/dbs/gen/tss/scmtni.py {input}/Homo_sapiens.GRCh37.74.TSS.5000.bed {output}
        """

rule gen_tss_promoters:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        img='workflow/envs/gretabench.sif',
        prom=rules.cre_promoters.output
    output: 'dbs/hg38/gen/tss/promoters.bed.gz'
    shell:
        """
        cp {input.prom} {output}
        """
