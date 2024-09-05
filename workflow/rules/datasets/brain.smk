SAMPLES = glob_wildcards('datasets/brain/{sample}.frags.tsv.gz').sample
print(SAMPLES)

rule download_brain:
    output:
        multi_smpl1=local('datasets/brain/multiome_original_C16007.h5'),
        frags_smpl1=local('datasets/brain/C16007.frags.tsv.gz'),
        multi_smpl2=local('datasets/brain/multiome_original_YC4345.h5'),
        frags_smpl2=local('datasets/brain/YC4345.frags.tsv.gz'),
        annot=local('datasets/brain/raw_annot.csv')

    params:
        multi_smpl1=config['datasets']['brain']['url']['multi_smpl1'],
        frags_smpl1=config['datasets']['brain']['url']['frags_smpl1'],
        multi_smpl2=config['datasets']['brain']['url']['multi_smpl2'],
        frags_smpl2=config['datasets']['brain']['url']['frags_smpl2'],
        annot=config['datasets']['brain']['url']['annot']

    shell:
        """
        wget '{params.frags_smpl1}' -O '{output.frags_smpl1}'
        wget '{params.multi_smpl1}' -O '{output.multi_smpl1}'
        wget '{params.frags_smpl2}' -O '{output.frags_smpl2}'
        wget '{params.multi_smpl2}' -O '{output.multi_smpl2}'
        wget '{params.annot}' -O '{output.annot}'

        """

rule prc_annot:
    input:
        raw_annot=local('datasets/brain/raw_annot.csv'),
        sample_dir=local(directory('datasets/brain/'))
    
    output:
        annot=local('datasets/brain/annot.csv')

    singularity:
        'workflow/envs/gretabench.sif'
        
    shell:
        """
        python workflow/scripts/datasets/brain/prc_annot.py \
        -a {input.raw_annot} \
        -b {output.annot} \
        -c {input.sample_dir} 
        """


rule callpeaks_brain:
    threads: 4
    input:
        frags=expand('datasets/brain/{sample}.frags.tsv.gz', sample=SAMPLES),
        annot='datasets/brain/annot.csv',
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        tmp=temp(directory(local('datasets/brain/tmp_peaks'))),
        peaks=local('datasets/brain/peaks.h5ad')
    resources:
        mem_mb=8000,
        runtime=2160,
    shell:
        """
        python workflow/scripts/datasets/callpeaks.py \
        -f {input.frags} \
        -a {input.annot} \
        -t {output.tmp} \
        -o {output.peaks}
        """