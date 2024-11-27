rule download_reprofibro:
    threads: 2
    singularity: 'workflow/envs/figr.sif'
    output:
        tar=temp(local('dts/reprofibro/RAW.tar')),
        barcodes=temp(local('dts/reprofibro/barcode_map.tsv.gz')),
        genes=temp(local('dts/reprofibro/genes.tsv.gz')),
        bars=temp(local(expand('dts/reprofibro/{sample}.barcodes.tsv.gz', sample=config['dts']['reprofibro']['samples']))),
        frags=expand('dts/reprofibro/{sample}.frags.tsv.gz', sample=config['dts']['reprofibro']['samples']),
        tbis=expand('dts/reprofibro/{sample}.frags.tsv.gz.tbi', sample=config['dts']['reprofibro']['samples']),
        mats=temp(local(expand('dts/reprofibro/{sample}.matrix.mtx.gz', sample=config['dts']['reprofibro']['samples']))),
        annot=temp(local('dts/reprofibro/annot.csv')),
    params:
        tar=config['dts']['reprofibro']['url']['tar'],
        barcodes=config['dts']['reprofibro']['url']['barcodes'],
        genes=config['dts']['reprofibro']['url']['genes'],
        annot=config['dts']['reprofibro']['url']['annot'],
    shell:
        """
        data_path=$(dirname {output.tar}) && \
        wget --no-verbose '{params.tar}' -O {output.tar} && \
        tar xvf {output.tar} -C $data_path && \
        for file in $data_path/GSM*_D*.frag.bed.gz; do
            base_name=$(basename "$file" .frag.bed.gz);
            new_file="${{base_name#*_}}.frags.tsv.gz";
            mv $file $data_path/$new_file
        done && \
        ls $data_path/*.frags.tsv.gz | xargs -n 1 -P {threads} bash workflow/scripts/dts/format_frags.sh && \
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
        python workflow/scripts/dts/reprofibro/prc_annot.py -a {output.annot} && \
        awk 'BEGIN {{FS=OFS=","}} NR==1 {{print $0; next}} {{print $2"_"$1,$2,$3}}' {output.annot} > {output.annot}.tmp && \
        mv {output.annot}.tmp {output.annot} && \
        wget --no-verbose '{params.barcodes}' -O {output.barcodes} && \
        wget --no-verbose '{params.genes}' -O {output.genes}
        """


rule callpeaks_reprofibro:
    threads: 32
    singularity: 'workflow/envs/gretabench.sif'
    input:
        frags=rules.download_reprofibro.output.frags,
        annot=rules.download_reprofibro.output.annot,
    output: peaks=temp(local('dts/reprofibro/peaks.h5ad'))
    resources: mem_mb=64000
    shell:
        """
        python workflow/scripts/dts/callpeaks.py \
        -f {input.frags} \
        -a {input.annot} \
        -t '/tmp/reprofibro/' \
        -n {threads} \
        -o {output.peaks}
        """


rule annotate_reprofibro:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        mats=rules.download_reprofibro.output.mats,
        bars=rules.download_reprofibro.output.bars,
        path_gsym=rules.download_reprofibro.output.genes,
        path_peaks=rules.callpeaks_reprofibro.output.peaks,
        path_annot=rules.download_reprofibro.output.annot,
        path_barmap=rules.download_reprofibro.output.barcodes,
        gid=rules.gen_gid_ensmbl.output,
    output: out='dts/reprofibro/annotated.h5mu'
    shell:
        """
        python workflow/scripts/dts/reprofibro/reprofibro.py \
        -a {input.mats} \
        -b {input.bars} \
        -e {input.path_gsym} \
        -f {input.path_peaks} \
        -g {input.path_annot} \
        -i {input.path_barmap} \
        -j {input.gid} \
        -l {output}
        """
