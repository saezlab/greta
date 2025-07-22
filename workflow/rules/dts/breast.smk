localrules: prcannot_breast

rule download_fragments_breast:
    threads: 7
    singularity: 'workflow/envs/figr.sif'
    output:
        tar=temp(local('dts/breast/fragments.tsv.bgz')),
        frag=expand('dts/breast/{sample}.frags.tsv.gz', sample=config['dts']['breast']['samples'])
    params:
        tar=config['dts']['breast']['url']['tar']
    shell:
        """
        data_path=$(dirname "{output.tar}")
        wget --no-verbose '{params.tar}' -O '{output.tar}'
        # Decompress and split by pool
        gunzip -c {output.tar} | \
        awk 'BEGIN {{OFS="\\t"}} {{
            split($4, arr, "_");           # arr[1]=barcode, arr[2]=pool id
            gsub("-1", "", arr[1]);        # remove trailing -1
            pool = "ID" arr[2];          
            new_barcode = pool "_" arr[1];
            print $1, $2, $3, new_barcode, $5 >> "dts/breast/" pool ".frags.tsv"
        }}'

        # Compress each file with bgzip
        for f in dts/breast/ID*.frags.tsv; do
            bgzip -f -c "$f" > "$f.gz"
            tabix -p bed "$f.gz"
            rm "$f"
            rm "$f.gz.tbi"
        done

        rm '{output.tar}'
        """

rule download_anndata_breast:
    threads: 1
    output:
        adata=temp(local('dts/breast/multiome_raw.h5ad')),
    params:
        adata=config['dts']['breast']['url']['anndata']
    shell:
        """
        wget --no-verbose '{params.adata}' -O '{output.adata}'
        """
        
rule prcannot_breast:
    threads: 1
    singularity:
        'workflow/envs/gretabench.sif'
    input: rules.download_anndata_breast.output.adata,
    output:
        annot=temp(local('dts/breast/annot.csv'))
    shell:
        """
        python workflow/scripts/dts/breast/breast_annot.py \
        -i {input} \
        -o {output.annot}
        """

rule callpeaks_breast:
    threads: 8
    singularity: 'workflow/envs/gretabench.sif'
    input:
        frags=rules.download_fragments_breast.output.frag,
        annot=rules.prcannot_breast.output.annot,
    output: peaks=temp(local('dts/breast/peaks.h5ad'))
    resources:
        mem_mb=128000,
        runtime=2160,
    shell:
        """
        python workflow/scripts/dts/callpeaks.py \
        -f {input.frags} \
        -a {input.annot} \
        -t $TMPDIR \
        -n {threads} \
        -o {output.peaks}
        """

rule annotate_breast:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        path_h5ad=rules.download_anndata_breast.output.adata,
        path_peaks=rules.callpeaks_breast.output.peaks,
        path_annot=rules.prcannot_breast.output.annot,
        gid=rules.gen_gid_ensmbl.output
    output: out='dts/breast/annotated.h5mu'
    shell:
        """
        python workflow/scripts/dts/breast/breast.py \
        -a {input.path_h5ad} \
        -b {input.path_peaks} \
        -c {input.path_annot} \
        -e {input.gid} \
        -f {output}
        """