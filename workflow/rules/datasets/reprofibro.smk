localrules: download_reprofibro


rule download_reprofibro:
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
        rename_files() {{
            local data_path=$1
            local sufix=$2
            local new_sufix=$3
            shopt -s nullglob
            local files=($data_path/GSM*_*.$sufix)
            shopt -u nullglob
            for file in ${{files[@]}}; do
                local new_name=$(basename $file | sed 's/GSM776338[0-9]*_//' | sed 's/\.$sufix$/.$new_sufix/')
                mv $file $data_path/$new_name
            done
        }}
        rename_files $data_path barcodes.tsv.gz barcodes.tsv.gz
        rename_files $data_path frag.tsv.gz frags.tsv.gz
        rename_files $data_path matrix.mtx.gz matrix.mtx.gz
        wget --no-verbose '{params.barcodes}' -O {output.barcodes}
        wget --no-verbose '{params.genes}' -O {output.genes}
        """


rule callpeaks_reprofibro:
    threads: 32
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
        -o {output.peaks}
        """


rule annotate_reprofibro:
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
        -a {input.path_matrix_d1m} \
        -b {input.path_barcodes_d1m} \
        -c {input.path_matrix_d2m} \
        -d {input.path_barcodes_d2m} \
        -e {input.path_gsym} \
        -f {input.path_peaks} \
        -g {input.path_annot} \
        -i {input.path_barmap} \
        -j {input.path_geneids} \
        -k {params.organism} \
        -l {output}
        """
