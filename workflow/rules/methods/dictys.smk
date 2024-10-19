localrules: download_motifs_dictys, download_gene_annotations_dictys, download_genomes_dictys, pre_dictys


rule download_genomes_dictys:
    output:
        d=directory('gdata/dictys'),
    shell:
        """
        python workflow/scripts/methods/dictys/download_genomes.py \
        -d {output.d}
        """


rule download_motifs_dictys:
    params:
        url_h="https://hocomoco11.autosome.org/final_bundle/hocomoco11/full/HUMAN/mono/HOCOMOCOv11_full_HUMAN_mono_homer_format_0.0001.motif"
    output:
        h="gdata/dictys/motifs.motif",
    shell:
        """
        wget -O {output.h} {params.url_h}
        """


rule download_gene_annotations_dictys:
    conda:
        '../../envs/dictys.yaml'
    params:
        url_gtf="http://ftp.ensembl.org/pub/release-107/gtf/homo_sapiens/Homo_sapiens.GRCh38.107.gtf.gz"
    output:
        gtf=temp(local('gdata/dictys/gene.gtf')),
        bed='gdata/dictys/gene.bed'
    shell:
        """
        wget -O {output.gtf}.gz {params.url_gtf} && \
        gunzip {output.gtf}.gz
        dictys_helper gene_gtf.sh {output.gtf} {output.bed}
        """


rule pre_dictys:
    input:
        mdata=rules.extract_case.output.mdata
    output:
        out='datasets/{dataset}/cases/{case}/runs/dictys.pre.h5mu'
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    shell:
        """
        cp {input.mdata} {output.out}
        """


rule p2g_dictys:
    threads: 1
    conda:
        '../../envs/dictys.yaml'
    input:
        pre=lambda wildcards: map_rules('pre', wildcards.pre),
        annotation=rules.download_gene_annotations_dictys.output.bed,
    output:
        d=temp(directory(local('datasets/{dataset}/cases/{case}/runs/{pre}_dictys_tmp/'))),
        out='datasets/{dataset}/cases/{case}/runs/{pre}.dictys.p2g.csv',
    params:
        ext=config['methods']['dictys']['ext']
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    shell:
        """
        mkdir {output.d}
        set +e
        timeout $(({resources.runtime}-20))m bash -c \
        'python workflow/scripts/methods/dictys/p2g.py \
        -d {input.pre} \
        -t {output.d} \
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
        frags=list_frags_files,
        motif=rules.download_motifs_dictys.output.h,
        genome=rules.download_genomes_dictys.output.d,
    output:
        d=temp(directory('datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.dictys_tmp')),
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
        annotation=rules.download_gene_annotations_dictys.output.bed,
    output:
        d=temp(directory('datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.dictys_tmp')),
        out='datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.dicyts.mdl.csv'
    params:
        d=config['methods']['dictys']['device'],
        ext=config['methods']['dictys']['ext']
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
        -p {output.out} \
        -t {threads} \
        -v {params.d} \
        -e {params.ext} \
        -g {input.annotation}
        if [ $? -eq 124 ]; then
            awk 'BEGIN {{ print "source,target,score,pval" }}' > {output.out}
        fi
        """
