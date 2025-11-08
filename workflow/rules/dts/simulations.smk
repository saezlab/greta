localrules: download_simulations

rule download_simulations:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        img='workflow/envs/gretabench.sif',
    output: directory('dts/simulations')
    params:
        url=config['url_sim'],
    shell:
        """
        path_sim=$(dirname {output})
        path_tmp=$path_sim/tmp.zip
        wget --no-verbose "{params.url}" -O $path_tmp
        python -c "import zipfile, os; os.makedirs('$path_sim', exist_ok=True); [zipfile.ZipFile('$path_tmp').extract(m, '$path_sim') for m in zipfile.ZipFile('$path_tmp').namelist()]"
        #unzip -o $path_tmp -d $path_sim
        rm -rf $path_tmp
        for f in {output}/seed_*; do
            base=$(basename "$f")
            path_dataset={output}/$base
            python workflow/scripts/dts/simulations/pre.py $path_dataset
        done
        """