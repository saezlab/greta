rule download_reprofibro:
    threads: 2
    singularity: 'workflow/envs/figr.sif'
    output:
        tar=temp(local('datasets/reprofibro/RAW.tar')),
        barcodes=temp(local('datasets/reprofibro/barcode_map.tsv.gz')),
        genes=temp(local('datasets/reprofibro/genes.tsv.gz')),
        bars=temp(local(expand('datasets/reprofibro/{sample}.barcodes.tsv.gz', sample=config['datasets']['reprofibro']['samples']))),
        frags=expand('datasets/reprofibro/{sample}.frags.tsv.gz', sample=config['datasets']['reprofibro']['samples']),
        tbis=expand('datasets/reprofibro/{sample}.frags.tsv.gz.tbi', sample=config['datasets']['reprofibro']['samples']),
        mats=temp(local(expand('datasets/reprofibro/{sample}.matrix.mtx.gz', sample=config['datasets']['reprofibro']['samples']))),
        annot=temp(local('datasets/reprofibro/annot.csv')),
    params:
        tar=config['datasets']['reprofibro']['url']['tar'],
        barcodes=config['datasets']['reprofibro']['url']['barcodes'],
        genes=config['datasets']['reprofibro']['url']['genes'],
        annot=config['datasets']['reprofibro']['url']['annot'],
    shell:
        """
        data_path=$(dirname {output.tar})
        wget --no-verbose '{params.tar}' -O {output.tar} && \
        tar xvf {output.tar} -C $data_path && \
        for file in $data_path/GSM*_D*.frag.bed.gz; do
            base_name=$(basename "$file" .frag.bed.gz);
            new_file="${{base_name#*_}}.frags.tsv.gz";
            mv $file $data_path/$new_file
        done && \
        ls $data_path/*.frags.tsv.gz | xargs -n 1 -P {threads} bash workflow/scripts/datasets/format_frags.sh && \
        for file in $data_path/GSM*_*.barcodes.tsv.gz; do
            base_name=$(basename "$file" .barcodes.tsv.gz);
            new_file="${{base_name#*_}}.barcodes.tsv.gz";
            mv $file $data_path/$new_file;
        done && \
        for file in $data_path/GSM*_*.matrix.mtx.gz; do
            base_name=$(basename "$file" .matrix.mtx.gz);
            new_file="${{base_name#*_}}.matrix.mtx.gz";
            mv $file $data_path/$new_file;
        done && \
        wget --no-verbose '{params.annot}' -O {output.annot} && \
        python workflow/scripts/datasets/reprofibro/prc_annot.py -a {output.annot} && \
        awk 'BEGIN {{FS=OFS=","}} NR==1 {{print $0; next}} {{print $2"_"$1,$2,$3}}' {output.annot} > {output.annot}.tmp && \
        mv {output.annot}.tmp {output.annot} && \
        wget --no-verbose '{params.barcodes}' -O {output.barcodes} && \
        wget --no-verbose '{params.genes}' -O {output.genes}
        """


rule callpeaks_reprofibro:
    threads: 32
    singularity:
        'workflow/envs/gretabench.sif'
    input:
        frags=rules.download_reprofibro.output.frags,
        annot=rules.download_reprofibro.output.annot,
    output:
        tmp=temp(directory(local('datasets/reprofibro/tmp'))),
        peaks=temp(local('datasets/reprofibro/peaks.h5ad'))
    resources:
        mem_mb=64000,
    shell:
        """
        python workflow/scripts/datasets/callpeaks.py \
        -f {input.frags} \
        -a {input.annot} \
        -t {output.tmp} \
        -n {threads} \
        -o {output.peaks}
        """


rule annotate_reprofibro:
    threads: 1
    singularity:
        'workflow/envs/gretabench.sif'
    input:
        mats=rules.download_reprofibro.output.mats,
        bars=rules.download_reprofibro.output.bars,
        path_gsym=rules.download_reprofibro.output.genes,
        path_peaks=rules.callpeaks_reprofibro.output.peaks,
        path_annot=rules.download_reprofibro.output.annot,
        path_barmap=rules.download_reprofibro.output.barcodes,
        path_geneids=rules.download_geneids.output.dr,
    output:
        out='datasets/reprofibro/annotated.h5mu'
    params:
        organism=config['datasets']['reprofibro']['organism']
    shell:
        """
        python workflow/scripts/datasets/reprofibro/reprofibro.py \
        -a {input.mats} \
        -b {input.bars} \
        -e {input.path_gsym} \
        -f {input.path_peaks} \
        -g {input.path_annot} \
        -i {input.path_barmap} \
        -j {input.path_geneids} \
        -k {params.organism} \
        -l {output}
        """
