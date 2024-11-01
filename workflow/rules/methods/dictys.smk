localrules: download_motifs_dictys, download_gene_annotations_dictys


rule download_genomes_dictys:
    conda: '../../envs/dictys.yaml'
    output:
        d=directory('gdata/dictys_genomes/'),
    shell:
        """
        dictys_helper genome_homer.sh hg38 {output.d}
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
    threads: 1
    conda:
        '../../envs/dictys.yaml'
    input:
        mdata=rules.extract_case.output.mdata
    output:
        tmp='datasets/{dataset}/cases/{case}/runs/dictys_pre_expr.tsv.gz',
        out='datasets/{dataset}/cases/{case}/runs/dictys.pre.h5mu',
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    shell:
        """
        python workflow/scripts/methods/dictys/pre.py \
        -m {input.mdata} \
        -t {output.tmp} \
        -o {output.out}
        """


rule p2g_dictys:
    threads: 1
    conda:
        '../../envs/dictys.yaml'
    container: None
    input:
        pre=lambda wildcards: map_rules('pre', wildcards.pre),
        annotation=rules.download_gene_annotations_dictys.output.bed,
    output:
        d=temp(directory(local('datasets/{dataset}/cases/{case}/runs/{pre}_dictys_tmp/'))),
        out='datasets/{dataset}/cases/{case}/runs/{pre}.dictys.p2g.csv',
    params:
        ext=config['methods']['dictys']['ext'] // 2,
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    shell:
        """
        mkdir -p {output.d}
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
    conda:
        '../../envs/dictys.yaml'
    container: None
    input:
        pre=lambda wildcards: map_rules('pre', wildcards.pre),
        p2g=lambda wildcards: map_rules('p2g', wildcards.p2g),
        frags=list_frags_files,
        motif=rules.download_motifs_dictys.output.h,
        genome=rules.download_genomes_dictys.output.d,
    output:
        d=temp(directory('datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.dictys_tmp')),
        out='datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.dictys.tfb.csv'
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    params:
        use_p2g=True,
    shell:
        """
        set +e
        timeout $(({resources.runtime}-20))m \
        bash workflow/scripts/methods/dictys/tfb.sh \
        --input_pre {input.pre} \
        --input_p2g {input.p2g} \
        --input_frags {input.frags} \
        --input_motif {input.motif} \
        --input_genome {input.genome} \
        --output_d {output.d} \
        --output_out {output.out} \
        --threads {threads} \
        --use_p2g {params.use_p2g} 
        if [ $? -eq 124 ]; then
            awk 'BEGIN {{ print "cre,tf,score" }}' > {output.out}
        fi
        """


rule mdl_dictys:
    threads: 32
    conda:
        '../../envs/dictys.yaml'
    container: None
    input:
        pre=lambda wildcards: map_rules('pre', wildcards.pre),
        p2g=lambda wildcards: map_rules('p2g', wildcards.p2g),
        tfb=lambda wildcards: map_rules('tfb', wildcards.tfb),
        annotation=rules.download_gene_annotations_dictys.output.bed,
    output:
        d=temp(directory('datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.dictys_tmp')),
        out='datasets/{dataset}/cases/{case}/runs/{pre}.{p2g}.{tfb}.dicyts.mdl.csv'
    params:
        ext=config['methods']['dictys']['ext'] // 2,
        n_p2g_links=config['methods']['dictys']['n_p2g_links'],
        device=config['method']['dictys']['device'],
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    shell:
        """
        set +e
        timeout $(({resources.runtime}-20))m \
        bash workflow/scripts/methods/dictys/mdl.sh \
        --output_d {output.d} \
        --pre_path {input.pre} \
        --p2g_path {input.p2g} \
        --tfb_path {input.tfb} \
        --annot {input.annotation} \
        --distance {params.ext} \
        --n_p2g_links {params.n_p2g_links} \
        --threads {threads} \
        --out_path {output.out} \
        --device {params.device}
        if [ $? -eq 124 ]; then
            awk 'BEGIN {{ print "source,target,score,pval" }}' > {output.out}
        fi
        """


rule mdl_o_dictys:
    threads: 32
    conda: '../../envs/dictys.yaml'
    container: None
    input:
        mdata=rules.extract_case.output.mdata,
        annotation=rules.download_gene_annotations_dictys.output.bed,
        frags=list_frags_files,
        motif=rules.download_motifs_dictys.output.h,
        genome=rules.download_genomes_dictys.output.d,
    output:
        tmp='datasets/{dataset}/cases/{case}/runs/o_dictys_pre_expr.tsv.gz',
        d=temp(directory(local('datasets/{dataset}/cases/{case}/runs/o_dictys_dictys_tmp/'))),
        pre=temp(local('datasets/{dataset}/cases/{case}/runs/o_dictys.pre.h5mu')),
        p2g=temp(local('datasets/{dataset}/cases/{case}/runs/o_dictys.o_dictys.p2g.csv')),
        tfb=temp(local('datasets/{dataset}/cases/{case}/runs/o_dictys.o_dictys.o_dictys.tfb.csv')),
        out='datasets/{dataset}/cases/{case}/runs/o_dictys.o_dictys.o_dictys.o_dictys.grn.csv',
    params:
        ext=config['methods']['dictys']['ext'] // 2,
        n_p2g_links=config['methods']['dictys']['n_p2g_links'],
        device=config['method']['dictys']['device'],
        use_p2g=False,
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    shell:
        """
        set +e
        timeout $(({resources.runtime}-20))m bash -c \
        'mkdir -p {output.d} && \
        python workflow/scripts/methods/dictys/pre.py \
        -m {input.mdata} \
        -t {output.tmp} \
        -o {output.pre} && \
        python workflow/scripts/methods/dictys/p2g.py \
        -d {output.pre} \
        -t {output.d} \
        -p {output.p2g} \
        -e {params.ext} \
        -g {input.annotation} && \
        bash workflow/scripts/methods/dictys/tfb.sh \
        --input_pre {output.pre} \
        --input_p2g {output.p2g} \
        --input_frags {input.frags} \
        --input_motif {input.motif} \
        --input_genome {input.genome} \
        --output_d {output.d} \
        --output_out {output.tfb} \
        --use_p2g {params.use_p2g} \
        --threads {threads} && \
        bash workflow/scripts/methods/dictys/mdl.sh \
        --output_d {output.d} \
        --pre_path {output.pre} \
        --p2g_path {output.p2g} \
        --tfb_path {output.tfb} \
        --annot {input.annotation} \
        --distance {params.ext} \
        --n_p2g_links {params.n_p2g_links} \
        --threads {threads} \
        --device {params.device} \
        --out_path {output.out}'
        if [ $? -eq 124 ]; then
            awk 'BEGIN {{ print "cre,gene,score" }}' > {output.out}
        fi
        """
