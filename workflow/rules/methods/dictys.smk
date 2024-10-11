localrules: download_motifs, download_gene_annotations, download_genomes


rule download_genomes:
    output:
        d=directory('gdata/genomes'),
    shell:
        """
        python workflow/scripts/methods/dictys/download_genomes.py \
        -d {output.d}
        """


rule download_motifs:
    params:
        url_h="https://hocomoco11.autosome.org/final_bundle/hocomoco11/full/HUMAN/mono/HOCOMOCOv11_full_HUMAN_mono_homer_format_0.0001.motif"
    output:
        h="aertslab/motifs.motif",
    shell:
        """
        wget -O {output.h} {params.url_h}
        """


rule download_gene_annotations:
    params:
        url_h="http://ftp.ensembl.org/pub/release-107/gtf/homo_sapiens/Homo_sapiens.GRCh38.107.gtf.gz"
    output:
        h="aertslab/genomes/gene.gtf",
    shell:
        """
        wget -O gene.gtf.gz {params.url_h} && \
        gunzip gene.gtf.gz
        """


rule pre_dictys:
    threads: 32
    input:
        mdata=rules.extract_case.output.mdata
    output:
        out='datasets/{dataset}/cases/{case}/runs/dictys.pre.h5mu'
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    shell:
        """
        python workflow/scripts/methods/dictys/pre.py \
        -i {input.mdata} \
        -o {output.out}
        """


rule p2g_dictys:
    threads: 32
    input:
        pre=lambda wildcards: map_rules('pre', wildcards.pre),
        wd='datasets/{dataset}/cases/{case}/runs/',
        annotation=rules.download_gene_annotations.output.h,
    output:
        out='datasets/{dataset}/cases/{case}/runs/{pre}.dictys.p2g.csv',
    params:
        ext=config['methods']['dictys']['ext']
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    shell:
        """
        set +e
        timeout $(({resources.runtime}-20))m bash -c \
        'python workflow/scripts/methods/dictys/p2g.py \
        -d {input.pre} \
        -w {input.wd} \
        -p {output.out} \
        -e {params.ext} \
        -g {input.annotation}'
        if [ $? -eq 124 ]; then
            awk 'BEGIN {{ print "cre,gene,score" }}' > {output.out}
        fi
        """


rule tfb_dictys:
    threads: 32
    input:
        pre=lambda wildcards: map_rules('pre', wildcards.pre),
        p2g=lambda wildcards: map_rules('p2g', wildcards.p2g),
        wd='datasets/{dataset}/cases/{case}/runs/',
        frags=list_frags_files,
        motif=rules.download_motifs.output.h,
        genome=rules.download_genomes.output.d,
    output:
        out='datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.dicyts.tfb.csv'
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    shell:
        """
        set +e
        timeout $(({resources.runtime}-20))m \
        python workflow/scripts/methods/dictys/tfb.py \
        -d {input.pre} \
        -k {input.p2g} \
        -w {input.wd} \
        -f {input.frags} \
        -p {output.out} \
        -t {threads} \
        -m {input.motif} \
        -g {input.genome}
        if [ $? -eq 124 ]; then
            awk 'BEGIN {{ print "cre,tf,score" }}' > {output.out}
        fi
        """


rule mdl_dictys:
    threads: 32
    input:
        pre=lambda wildcards: map_rules('pre', wildcards.pre),
        p2g=lambda wildcards: map_rules('p2g', wildcards.p2g),
        tfb=lambda wildcards: map_rules('tfb', wildcards.tfb),
        wd='datasets/{dataset}/cases/{case}/runs/',
    output:
        out='datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.dicyts.mdl.csv'
    params:
        d=config['methods']['dictys']['device'],
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    shell:
        """
        set +e
        timeout $(({resources.runtime}-20))m \
        python workflow/scripts/methods/dictys/mdl.py \
        -d {input.pre} \
        -g {input.p2g} \
        -b {input.tfb} \
        -w {input.wd} \
        -p {output.out} \
        -t {threads} \
        -v {params.d}
        if [ $? -eq 124 ]; then
            awk 'BEGIN {{ print "source,target,score,pval" }}' > {output.out}
        fi
        """

