localrules: dwn_image

rule dwn_image:
    threads: 1
    singularity: None
    output: 'workflow/envs/{name_img}.sif'
    resources:
        mem_mb=8000,
        runtime=config['max_mins_per_step'],
    params: id=config['zenodo_id']
    shell:
        """
        wget "https://zenodo.org/records/{params.id}/files/{wildcards.name_img}.sif?download=1" -O {output}
        """
