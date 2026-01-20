localrules: prcannot_kidney, annotate_kidney

rule download_kidney:
    threads: 12
    singularity: 'workflow/envs/figr.sif'
    input: 'workflow/envs/figr.sif'
    output:
        h5se=temp(local('dts/hg38/kidney/data.h5se')),
        frag=temp(local(expand('dts/hg38/kidney/{sample}.frags.tsv.gz', sample=config['dts']['kidney']['samples'])))
    params:
        url_h5se=config['dts']['kidney']['url']['h5se'],
        url_tar=config['dts']['kidney']['url']['tar'],
        samples=config['dts']['kidney']['samples'],
    shell:
        """
        path_kidney=$(dirname {output.h5se})
        wget --no-verbose "{params.url_h5se}" -O {output.h5se}
        for SAMPLE in {params.samples}; do
            echo $SAMPLE
            path_zip=$(echo $path_kidney/${{SAMPLE}}.zip)
            wget --no-verbose "{params.url_tar}/KPMP_${{SAMPLE}}_10X_Dual.zip?download=1" -O $path_zip
            python -c "import zipfile; zipfile.ZipFile('$path_zip').extractall('$path_kidney')" && rm $path_zip
            path_sample=$(echo "$path_kidney/KPMP_${{SAMPLE}}_10X_Dual")
            mv "${{path_sample}}/atac_fragments.tsv.gz" "$path_kidney/${{SAMPLE}}.frags.tsv.gz"
            rm -r $path_sample
        done
        ls $path_kidney/*.frags.tsv.gz | xargs -n 1 -P {threads} bash workflow/scripts/dts/format_frags.sh
        """


rule prcannot_kidney:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        h5se=rules.download_kidney.output.h5se
    output:
        annot=temp(local('dts/hg38/kidney/annot.csv')),
        rna=temp(local('dts/hg38/kidney/rna.h5ad'))
    params:
        samples=config['dts']['kidney']['samples'],
    shell:
        """
        python workflow/scripts/dts/kidney/annot.py \
        -a {input.h5se} \
        -b {output.annot} \
        -c {output.rna} \
        -d {params.samples}
        """


rule callpeaks_kidney:
    threads: 12
    singularity: 'workflow/envs/gretabench.sif'
    input:
        frags=rules.download_kidney.output.frag,
        annot=rules.prcannot_kidney.output.annot,
    output: peaks=temp(local('dts/hg38/kidney/peaks.h5ad'))
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

rule annotate_kidney:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        gex=rules.prcannot_kidney.output.rna,
        peaks=rules.callpeaks_kidney.output.peaks,
        annot=rules.prcannot_kidney.output.annot,
        gid=rules.gen_gid_ensmbl.output.hg38,
    output: out='dts/hg38/kidney/annotated.h5mu'
    shell:
        """
        python workflow/scripts/dts/kidney/kidney.py \
        -a {input.gex} \
        -b {input.peaks} \
        -c {input.annot} \
        -e {input.gid} \
        -f {output.out}
        """
