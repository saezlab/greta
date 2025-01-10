## Adding databases

Databases are sorted into different categories:
- c2g: CRE to gene databases (e.g. eQTL studies)
- cre: CRE annotation (e.g. ENCODE)
- gen: general genome annotations (e.g. lambert, ENSMBL)
- gst: Gene set databases (e.g. REACTOME)
- ont: Ontology databases
- prt: TF perturbation databases (e.g. KnockTF)
- tfb: TF binding databases (e.g. ChIP-Atlas)
- tfm: TF marker databases (e.g. TF-Marker)
- tfp: TF-TF interaction databases (e.g. IntAct)
- tss: TSS databases (e.g. ENSMBL)

A url where to download the database should be provided in the `config/config.yaml` file.
```
# Databases
dbs:
    hg38:
        gen:
            ...
        prt:
            ...
        gst:
            ...
            newdatabase: 'https:// ...'
    mm10:
        ...
```
Note that databases are divided by organism.

Rules for each of these categories can be found in `workflow/rules/dbs/`. New databases should be added to their corresponding rule file.
If a database does not fit any of these categories, a new rule file can be created.

Here is an example of a rule for CRE:
```
rule cre_encode:
    threads: 1
    output: 'dbs/hg38/cre/encode/encode.bed'
    params:
        url=config['dbs']['hg38']['cre']['encode']
    shell:
        """
        ...
        """
```

Rules should follow this naming convetion: `{dbtype_dbname}`, in this case `cre_encode`.
The output should be stored using this path format: `dbs/{organism}/{dbtype}/{dbname}/{dbname}.bed`
When possible use `.bed` format, else `csv`.
