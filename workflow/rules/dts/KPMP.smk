localrules: prcannot_KPMP

rule download_fragments_KPMP:
    threads: 7
    singularity: 'workflow/envs/figr.sif'
    output:
        tar=temp(local('dts/KPMP/8029990.zip')),
        frag=expand('dts/KPMP/{sample}.frags.tsv.gz', sample=config['dts']['KPMP']['samples']),
        adata=temp(local('dts/KPMP/multiome_raw.h5Seurat'))
    params:
        tar=config['dts']['KPMP']['url']['tar']
    shell:
        """
        wget --no-verbose '{params.tar}' -O '{output.tar}'
        python3 -c "import zipfile; zipfile.ZipFile('{output.tar}').extractall('dts/KPMP/')"
        echo "unzipped multiome file"
        mv dts/KPMP/raw_snATAC_and_snRNA_human_kidney.h5Seurat '{output.adata}'
        echo "renamed adata"
        rm '{output.tar}'

        workflow/scripts/dts/KPMP/KPMP_process_frags.sh
        """

        
rule convert_anndata_KPMP:
    threads: 1
    input: rules.download_fragments_KPMP.output.adata
    output:
        adata=temp(local('dts/KPMP/multiome_raw.h5ad')),
    conda:
        '../../../workflow/envs/seuratdisk.yaml'
    shell:
        """
        Rscript -e "library(SeuratDisk); Convert('{input}', dest = 'h5ad')"
        """

rule prcannot_KPMP:
    threads: 1
    singularity:
        'workflow/envs/gretabench.sif'
    input: rules.convert_anndata_KPMP.output.adata,
    output:
        annot=temp(local('dts/KPMP/annot.csv'))
    shell:
        """
        python workflow/scripts/dts/KPMP/KPMP_annot.py \
        -i {input} \
        -o {output.annot}
        """

rule callpeaks_KPMP:
    threads: 8
    singularity: 'workflow/envs/gretabench.sif'
    input:
        frags=rules.download_fragments_KPMP.output.frag,
        annot=rules.prcannot_KPMP.output.annot,
    output: peaks=temp(local('dts/KPMP/peaks.h5ad'))
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

rule annotate_KPMP:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        path_h5ad=rules.convert_anndata_KPMP.output.adata,
        path_peaks=rules.callpeaks_KPMP.output.peaks,
        path_annot=rules.prcannot_KPMP.output.annot,
        gid=rules.gen_gid_ensmbl.output
    output: out='dts/KPMP/annotated.h5mu'
    shell:
        """
        python workflow/scripts/dts/KPMP/KPMP.py \
        -a {input.path_h5ad} \
        -b {input.path_peaks} \
        -c {input.path_annot} \
        -e {input.gid} \
        -f {output}
        """