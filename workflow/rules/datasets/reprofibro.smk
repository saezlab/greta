rule download_reprofibro:
    threads: 1
    output:
        tar=temp(local('datasets/reprofibro/RAW.tar')),
        barcodes=temp(local('datasets/reprofibro/barcode_map.tsv.gz')),
        genes=temp(local('datasets/reprofibro/genes.tsv.gz')),
        bars=temp(local(expand('datasets/reprofibro/{sample}.barcodes.tsv.gz', sample=config['datasets']['reprofibro']['samples']))),
        frags=expand('datasets/reprofibro/{sample}.frags.tsv.gz', sample=config['datasets']['reprofibro']['samples']),
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
        wget --no-verbose '{params.annot}' -O {output.annot}
        python workflow/scripts/datasets/reprofibro/prc_annot.py -a {output.annot}
        wget --no-verbose '{params.tar}' -O {output.tar}
        tar xvf {output.tar} -C $data_path
        cd $data_path
        for file in GSM*_*.frag.bed.gz; do
            new_file=$(echo $file | sed -E 's/GSM[0-9]+_([A-Za-z0-9]+)\.frag\.bed\.gz/\\1.frags.tsv.gz/')
            mv $file $new_file
            bash workflow/scripts/datasets/format_frags.sh $new_file
        done
        for file in GSM*_*.barcodes.tsv.gz; do
            new_file=$(echo $file | sed -E 's/GSM[0-9]+_([A-Za-z0-9]+)\.barcodes\.tsv\.gz/\\1.barcodes.tsv.gz/')
            mv $file $new_file
        done
        for file in GSM*_*.matrix.mtx.gz; do
            new_file=$(echo $file | sed -E 's/GSM[0-9]+_([A-Za-z0-9]+)\.matrix\.mtx\.gz/\\1.matrix.mtx.gz/')
            mv $file $new_file
        done
        cd ../../
        wget --no-verbose '{params.barcodes}' -O {output.barcodes}
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
