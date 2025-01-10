## Adding datasets

Datasets should contain both gene expression and chromatin accessibility (can be paired or unpaired). Namely, gene count matrices and fragment files are requiered. The annotation of barcodes into cell types and batches is also requiered.

Assuming the name of the dataset is newdataset, a rule file must be added in `workflow/rules/dts/newdataset.smk`.
For each rule, scripts running the specific steps should be placed in `workflow/scripts/dts/newdataset/example.py`.
For reference, other `.smk` files and associated scripts can be consulted for the already implemented datasets.
To ensure reproducibility, use the provided singularity images `workflow/envs/gretabench.sif` and `workflow/envs/figr.sif` (needed to format fragment files).
The following steps need to be included:

### 1- Adding metadata
In `config/config.yaml` the following information should be added:
```
# Datasets
dts:
    pbmc10k:
        ...
    newdataset:
        organism: 'hg38'  # Needed to process peaks and databases when evaluating
        url:
            data: 'https:// ...'  # This can be skipped if data does not need to be downloaded
        samples:
            ['sample_a', 'sample_b', ...]  # Sample (batch) information, this is used later on for integration
        cases:
            all:
                celltypes: 'all'  # Which cell types to use, here all, else pass a list of cell types to subset
                n_hvg: 16384  # Number of genes to consider
                n_hvr: 65536  # Number of peaks to consider
```

### 2- Data download and formatting of fragment files
In this rule, the download of the data is defined (can be a local path if already downloaded), and a formating of the fragment files is performed.
```
rule download_newdataset:
    threads: 8
    singularity: 'workflow/envs/figr.sif'  # Neded for format_frags.sh
    output:
        frag=expand('dts/heartatlas/{sample}.frags.tsv.gz', sample=config['dts']['newdataset']['samples']),  # Generates one fragment file per sample
        tbis=expand('dts/heartatlas/{sample}.frags.tsv.gz.tbi', sample=config['dts']['newdataset']['samples'])  # Generates one indexed fragment file per sample
    params:
        url=config['dts']['newdataset']['url']['data']  # Either pass url to download, or instead add input: path/to/file
    shell:
        """
        data_path=$(dirname "{output.frag[0]}")
        wget --no-verbose '{params.url}' ...  # Download data
        ls $data_path/*.frags.tsv.gz | xargs -n 1 -P {threads} bash workflow/scripts/dts/format_frags.sh  # Needed to format fragment files
```

### 3- Process barcode annotation
Barcode annotation should be extracted from the original data and processed. How to do this will depend on how the raw data for the new dataset was structured.
```
rule prcannot_newdataset:
    threads: 1
    input:
        ...
    output:
        annot=temp(local('dts/newdataset/annot.csv'))  # Where to store the barcode metadata
    shell:
        """
        ...
        """
```
The barcode metadata should have this format:
```
barcode,celltype,batch
sample1_AAACAGCCAACCGCCA,Stem cells,sample1
sample2_AAACAGCCACATTGCA,Gonadotropes,sample2
...
```
The barcode format should be `{sample_id}_{barcode}`
Note: There should not be any trailing -1 like in `AAACAGCCAACCGCCA-1` added in some single-cell processing frameworks.

### 4- Call chromatin accessibility peaks
This rule calls peaks from fragment files and the barcode annotation:
```
rule callpeaks_heartatlas:
    threads: 32
    singularity: 'workflow/envs/gretabench.sif'
    input:
        frags=rules.download_newdataset.output.frag,
        annot=rules.prcannot_newdataset.output.annot,
    output: peaks=temp(local('dts/newdataset/peaks.h5ad'))
    resources:
        mem_mb=64000,  # 64GBs, change as necessary
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
```

### 5- Annotate and clean final dataset
The obtained gene and peak count matrices have to be annotated by their barcode metadata and merged into a single [`MuData` object](https://mudata.readthedocs.io/en/latest/):
```
rule annotate_newdataset:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        ...
        path_peaks=rules.callpeaks_heartatlas.output.peaks,
        path_annot=rules.prcannot_heartatlas.output.annot,
        gid=lambda w: "dbs/{config['dts']['newdataset']['organism']}/gen/gid/ensembl.csv",
    output: out='dts/newdataset/annotated.h5mu'
    shell:
        """
        python workflow/scripts/dts/newdataset/newdataset.py ...
        """
```
The final mdata object should have this format:
```
MuData object with n_obs × n_vars = XXX × XXX
  obs:	'celltype', 'batch'
  2 modalities
    rna:	XXX x XXX
    atac:	XXX x XXX
```
Where inside each `AnnData` object (`rna` and `atac`), their `.X` attribute contains the raw integer counts (example: [0, 6, 2, 3, ...]). Normalization, integration and other qc steps are taken care by the pipeline automatically across downstream rules.



