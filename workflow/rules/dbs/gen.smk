localrules: gen_tfs_lambert, gen_tfs_scenic
localrules: gen_gid_ensmbl, gen_pid_uniprot, gen_genome_celloracle, gen_genome_dictys
localrules: gen_genome_scenicplus
localrules: gen_ann_dictys, gen_ann_pando
localrules: gen_motif_granie, gen_motif_dictys, gen_motif_scenic_rnk, gen_motif_scenic, gen_motif_scmtni
localrules: gen_motif_scenicplus, gen_genome_crema, gen_motif_crema, gen_genome_inferelator, download_liftover_chains


rule install_dictys:
    conda: '../../envs/dictys.yaml'
    output: 'workflow/envs/dictys.txt'
    shell:
        """
        pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118 && \
        touch {output}
        """


rule gen_tfs_lambert:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: 'workflow/envs/gretabench.sif'
    output: 'dbs/hg38/gen/tfs/lambert.csv'
    params: url=config['dbs']['hg38']['gen']['lambert']
    shell: "wget --no-check-certificate --no-verbose '{params.url}' -O {output}"


rule gen_tfs_scenic:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: 'workflow/envs/gretabench.sif'
    output: 'dbs/hg38/gen/tfs/scenic.csv'
    params:
        url=config['dbs']['hg38']['gen']['scenic']
    shell:
        "wget --no-verbose '{params.url}' -O {output}"


rule gen_gid_ensmbl:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: 'workflow/envs/gretabench.sif'
    output: expand('dbs/{org}/gen/gid/ensembl.csv', org=orgms)
    shell: "Rscript workflow/scripts/dbs/gen/gid/ensmbl.R {output}"


rule gen_pid_uniprot:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: 'workflow/envs/gretabench.sif'
    output: expand('dbs/{org}/gen/pid/uniprot.csv', org=orgms)
    shell: "Rscript workflow/scripts/dbs/gen/pid/uniprot.R {output}"


rule gen_genome_celloracle:
    threads: 4
    singularity: 'workflow/envs/celloracle.sif'
    input: 'workflow/envs/celloracle.sif'
    output: directory('dbs/hg38/gen/genome/celloracle/')
    shell:
        """
        python workflow/scripts/dbs/gen/genome/celloracle.py -o {output} &&
        mv {output}/hg38/* {output} && rm -r {output}/hg38/
        """


rule gen_genome_crema:
    threads: 1
    singularity: 'workflow/envs/crema.sif'
    input: 'workflow/envs/crema.sif'
    output: 'dbs/hg38/gen/genome/crema/genome.bed'
    shell:
        """
        wget 'https://raw.githubusercontent.com/zidongzh/CREMA/refs/heads/main/crema_resource/genomes/human.hg38.genome.bed' -O {output}
        """


rule gen_genome_dictys:
    threads: 4
    conda: '../../envs/dictys.yaml'
    input: rules.install_dictys.output
    output: directory('dbs/hg38/gen/genome/dictys/')
    shell:
        """
        dictys_helper genome_homer.sh hg38 {output}
        """


rule gen_genome_scenicplus:
    threads: 4
    singularity: 'workflow/envs/scenicplus.sif'
    input: 'workflow/envs/scenicplus.sif'
    output:
        ann='dbs/hg38/gen/genome/scenicplus/annotation.tsv',
        csz='dbs/hg38/gen/genome/scenicplus/chromsizes.tsv',
        tss='dbs/hg38/gen/genome/scenicplus/tss.tsv',
    shell:
        """
        scenicplus prepare_data download_genome_annotations \
        --species hsapiens \
        --genome_annotation_out_fname {output.ann} \
        --chromsizes_out_fname {output.csz}
        pycistopic tss get_tss \
        --output {output.tss} \
        --name "hsapiens_gene_ensembl" \
        --to-chrom-source ucsc \
        --ucsc hg38 \
        --no-cache
        """

rule gen_genome_inferelator:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: 'workflow/envs/gretabench.sif'
    output:
        gtf='dbs/hg38/gen/genome/inferelator/annotation.gtf.gz',
        tfa=temp(local('dbs/hg38/gen/genome/inferelator/tmp_genome.fa')),
        fa='dbs/hg38/gen/genome/inferelator/genome.fa.gz',
    shell:
        """
        wget --no-verbose 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz' -O {output.gtf}
        wget --no-verbose 'http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz' -O {output.tfa}.gz
        gunzip {output.tfa}.gz
        bgzip {output.tfa}
        cp {output.tfa}.gz {output.fa}
        mv {output.tfa}.gz {output.tfa}
        """


rule gen_motif_crema:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: 'workflow/envs/gretabench.sif'
    output: 'dbs/hg38/gen/motif/crema/jaspar.rds'
    shell:
        """
        wget 'https://github.com/zidongzh/CREMA/raw/refs/heads/main/crema_resource/tf_motifs/data_rds/motifs_pwmlist_human_jaspar.rds' -O {output}
        """



rule gen_motif_granie:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: 'workflow/envs/gretabench.sif'
    output: directory('dbs/hg38/gen/motif/granie/')
    shell:
        """
        wget --no-verbose 'https://s3.embl.de/zaugg-web/GRaNIE/TFBS/hg38/PWMScan_HOCOMOCOv12_H12INVIVO.tar.gz' -O {output}.tar.gz && \
        mkdir {output} && \
        tar -xvf {output}.tar.gz -C {output} && \
        find {output} -type f -exec mv {{}} {output} \; && \
        find {output} -type d -empty -delete && \
        rm {output}.tar.gz
        """


rule gen_motif_dictys:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: 'workflow/envs/gretabench.sif'
    params: url="https://hocomoco11.autosome.org/final_bundle/hocomoco11/full/HUMAN/mono/HOCOMOCOv11_full_HUMAN_mono_homer_format_0.0001.motif"
    output: 'dbs/hg38/gen/motif/dictys/dictys.motif'
    shell:
        """
        wget --no-verbose {params.url} -O {output}
        """


rule gen_motif_scenic_rnk:
    threads: 1
    singularity: 'workflow/envs/scenicplus.sif'
    input: 'workflow/envs/scenicplus.sif'
    output:
        sml='dbs/hg38/gen/motif/scenic/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather',
        big='dbs/hg38/gen/motif/scenic/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather'
    params:
        sml='https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather',
        big='https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather'
    shell:
        """
        wget --no-verbose {params.sml} -O {output.sml}
        wget --no-verbose {params.big} -O {output.big}
        """


rule gen_motif_scenic:
    threads: 1
    singularity: 'workflow/envs/scenicplus.sif'
    input: 'workflow/envs/scenicplus.sif'
    params:
        url="https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
    output: "dbs/hg38/gen/motif/scenic/nr.hgnc-m0.001-o0.0.tbl"
    shell:
        """
        wget --no-verbose {params.url} -O {output}
        """


rule gen_motif_scenicplus:
    threads: 1
    singularity: 'workflow/envs/scenicplus.sif'
    input: 'workflow/envs/scenicplus.sif'
    output:
        human_rankings="dbs/hg38/gen/motif/scenicplus/human_motif_SCREEN.regions_vs_motifs.rankings.feather",
        human_scores="dbs/hg38/gen/motif/scenicplus/human_motif_SCREEN.regions_vs_motifs.scores.feather",
        mouse_rankings="dbs/mm10/gen/motif/scenicplus/mouse_motif_SCREEN.regions_vs_motifs.rankings.feather",
        mouse_scores="dbs/mm10/gen/motif/scenicplus/mouse_motif_SCREEN.regions_vs_motifs.scores.feather",
        human_annot="dbs/hg38/gen/motif/scenicplus/motifs-v10nr_clust/nr.hgnc-m0.001-o0.0.tbl",
        mouse_annot="dbs/mm10/gen/motif/scenicplus/motifs-v10nr_clust/nr.mgi-m0.001-o0.0.tbl",
    params:
        human_rankings_url="https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/screen/mc_v10_clust/region_based/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather",
        human_scores_url="https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/screen/mc_v10_clust/region_based/hg38_screen_v10_clust.regions_vs_motifs.scores.feather",
        mouse_rankings_url="https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/screen/mc_v10_clust/region_based/mm10_screen_v10_clust.regions_vs_motifs.rankings.feather",
        mouse_scores_url="https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/screen/mc_v10_clust/region_based/mm10_screen_v10_clust.regions_vs_motifs.scores.feather",
        human_annot_url="https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl",
        mouse_annot_url="https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl",
    shell:
        """
        wget --no-verbose -O {output.human_rankings} {params.human_rankings_url}
        wget --no-verbose -O {output.human_scores} {params.human_scores_url}
        wget --no-verbose -O {output.mouse_rankings} {params.mouse_rankings_url}
        wget --no-verbose -O {output.mouse_scores} {params.mouse_scores_url}
        wget --no-verbose -O {output.human_annot} {params.human_annot_url}
        wget --no-verbose -O {output.mouse_annot} {params.mouse_annot_url}
        """


rule gen_motif_scmtni:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: 'workflow/envs/gretabench.sif'
    output: 
        motifs_dir = directory('dbs/hg38/gen/motif/scmtni/motif_files'),
        promoters_dir = directory('dbs/hg38/gen/motif/scmtni/promoter_files')
    params: 
        url="https://zenodo.org/records/8323399/files/RawMotifFiles.tar.gz?download=1"
    shell:
        """
        path_raw=$(dirname {output.motifs_dir})
        wget --no-verbose {params.url} -O $path_raw.tar.gz && \
        mkdir -p $path_raw && \
        tar -xvf $path_raw.tar.gz -C $path_raw && \
        path_raw=$path_raw/RawMotifFiles && \
        mkdir -p {output.motifs_dir} {output.promoters_dir} && \
        mv $path_raw/human_all_motifs_sorted_clean.txt {output.motifs_dir}/  && \
        mv $path_raw/Homo_sapiens.GRCh37.74.TSS.5000.bed {output.promoters_dir}/ && \
        rm $path_raw.tar.gz && \
        rm -rf $path_raw
        """


rule download_liftover_chains:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: 'workflow/envs/gretabench.sif'
    output:
        hg38ToHg19 = "dbs/chain_files/hg38ToHg19.over.chain.gz"
    params:
        url_hg38ToHg19 = "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz"
    shell:
        """
        mkdir -p dbs/chain_files
        wget --no-verbose {params.url_hg38ToHg19} -O {output.hg38ToHg19}
        """



rule gen_motif_inferelator:
    threads: 4
    singularity: 'workflow/envs/inferelator.sif'
    output:
        zip=temp(local('dbs/hg38/gen/motif/inferelator/pwm.zip')),
        dir=directory('dbs/hg38/gen/motif/inferelator/pwm'),
        meme='dbs/hg38/gen/motif/inferelator/cisbp.meme'
    shell:
        """
        wget --no-verbose 'https://cisbp.ccbr.utoronto.ca/tmp/Homo_sapiens_2025_08_22_6:07_pm.zip' -O {output.zip}
        unzip {output.zip} -d {output.dir}
        python -m inferelator_prior.pwm_to_meme --motif {output.dir}/pwms_all_motifs/* --info {output.dir}/TF_Information_all_motifs_plus.txt --out {output.meme}
        """


rule gen_ann_dictys:
    threads: 4
    conda: '{home_path}/miniforge3/envs/dictys'.format(home_path=home_path)
    params:
        url="http://ftp.ensembl.org/pub/release-107/gtf/homo_sapiens/Homo_sapiens.GRCh38.107.gtf.gz"
    input: 'workflow/envs/gretabench.sif'
    output: 'dbs/hg38/gen/ann/dictys/ann.bed'
    shell:
        """
        wget --no-verbose {params.url} -O {output}.gtf.gz && \
        gunzip {output}.gtf.gz
        dictys_helper gene_gtf.sh {output}.gtf {output} && \
        rm {output}.gtf
        """


rule gen_ann_pando:
    threads: 1
    singularity: 'workflow/envs/pando.sif'
    input: 'workflow/envs/pando.sif'
    output: 'dbs/hg38/gen/ann/pando/ann.csv'
    shell:
        """
        Rscript -e 'library(EnsDb.Hsapiens.v86); \
        gene.ranges_hg <- Signac::GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86); \
        write.csv(gene.ranges_hg, "{output}", row.names=FALSE);'
        """

