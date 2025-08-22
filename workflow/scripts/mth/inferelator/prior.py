import os
import pandas as pd
import numpy as np
import mudata as mu
from inferelator_prior import network_from_motifs as nfm
import sys
from inferelator_prior.motifs.motif_scan import MotifScan
from inferelator_prior.processor.prior import MotifScorer
from tqdm import tqdm
import pyranges as pr


# Read
path_data = sys.argv[1]
path_fa = sys.argv[2]
path_gtf = sys.argv[3]
path_motif = sys.argv[4]
window_size = int(sys.argv[5])
path_out = sys.argv[6]

# Read
mdata = mu.read(path_data)
genes = mdata.mod['rna'].var_names.astype('U')
peaks = mdata.mod['atac'].var_names.astype('U')

# Format
peaks = pd.DataFrame(peaks, columns=['cre'])
peaks[['Chromosome', 'Start', 'End']] = peaks['cre'].str.split('-', n=2, expand=True)
peaks = peaks[['Chromosome', 'Start', 'End']]

# Write
celltypes = np.unique(mdata.obs['celltype'].values)
path_adatas = os.path.join(path_out, 'adatas')
for i, celltype in enumerate(celltypes):
    tmp = mdata.mod['rna'][mdata.obs['celltype'] == celltype]
    os.makedirs(path_adatas, exist_ok=True)
    celltype = celltype.strip().replace(' ', '_')
    tmp.write(os.path.join(path_adatas, f'task_{celltype}.h5ad'))
pd.DataFrame(genes.values.reshape(-1, 1)).to_csv(os.path.join(path_out, 'genes.tsv'), header=False, index=False, sep='\t')
path_peaks = os.path.join(path_out, 'peaks.tsv')
peaks.to_csv(path_peaks, header=False, index=False, sep='\t')

# Scan
scanner_type='fimo'
MotifScan.set_type(scanner_type)
print('Loading motifs')
motifs, motif_information = nfm.load_and_process_motifs(
    path_motif,
    motif_format='meme',
    truncate_prob=0.35,
    regulator_constraint_list=genes,
    fuzzy=False,
    shuffle=None,
)
motif_information = motif_information[motif_information[nfm.MOTIF_NAME_COL].isin(genes)]
gene_constraint_list = list(motif_information['Motif_Name'].unique())
print('Motifs loaded')

print("Loading genes from file ({f})".format(f=path_gtf))
fasta_gene_len = nfm.get_fasta_lengths(path_fa)
genes = nfm.load_gtf_to_dataframe(
    path_gtf,
    fasta_record_lengths=fasta_gene_len
)
print(f"{genes.shape[0]} genes loaded")

genes = nfm.select_genes(genes, gene_constraint_list)
genes = nfm.open_window(
    genes,
    window_size=window_size,
    use_tss=True,
    fasta_record_lengths=fasta_gene_len,
    constrain_to_intergenic=None,
)
print(
    f"Promoter (n={genes.shape[0]}) regions defined with window {window_size} "
    f"around TSS"
)

pr_genes = genes[[nfm.GTF_CHROMOSOME, nfm.SEQ_START, nfm.SEQ_STOP, 'gene_name']].rename(columns={'seqname': 'Chromosome', 'start': 'Start', 'end': 'End'})
pr_genes = pr.PyRanges(pr_genes)
pr_peaks = pr.PyRanges(peaks)
p2g = pr_genes.intersect(pr_peaks).df.rename(columns={'gene_name': 'gene'})
p2g['cre'] = p2g['Chromosome'].astype(str) + '-' + p2g['Start'].astype(str) + '-' + p2g['End'].astype(str)
p2g['score'] = 1
p2g = p2g[['cre', 'gene', 'score']]
p2g.to_csv(os.path.join(path_out, 'p2g.csv'), index=False)

gene_locs = genes.loc[
    :,
    [nfm.GTF_CHROMOSOME, nfm.SEQ_START, nfm.SEQ_STOP, nfm.GTF_STRAND]
].copy()
gene_locs[
    [nfm.SEQ_START, nfm.SEQ_STOP]
] = gene_locs[[nfm.SEQ_START, nfm.SEQ_STOP]].astype(int)
extract_fasta = MotifScan.scanner.extract_genome(
    path_fa,
    constraint_bed_file=path_peaks,
    promoter_bed=gene_locs,
)

def network_scan_build_single_tf(
    tf_mi_df,
    motif_ic=None,
):
    tf_motifs = tf_mi_df[nfm.MOTIF_OBJ_COL].tolist()
    motif_peaks = MotifScan.scanner(
        motifs=tf_motifs,
        num_workers=1
    ).scan(
        None,
        extracted_genome=extract_fasta,
        min_ic=motif_ic,
        threshold="1e-4",
    )
    MotifScorer.set_information_criteria(
        min_binding_ic=motif_ic,
        max_dist=100
    )
    if motif_peaks is not None:
        ra_ma, pr_da = nfm.summarize_target_per_regulator(
            genes,
            motif_peaks,
            tf_mi_df,
            num_workers=1,
            silent=True,
            by_chromosome=True,
        )
    else:
        ra_ma = pd.DataFrame(
            0.0,
            index=genes[nfm.GTF_GENENAME],
            columns=tf_mi_df[nfm.MOTIF_NAME_COL].unique().tolist(),
        )
        pr_da = None
    net_results = nfm.network_build(
        ra_ma,
        pr_da,
        num_cores=1,
        output_prefix=None,
        silent=True
    )
    return net_results[0]

prior_matrix = []
for _, tf_mi_df in tqdm(motif_information.groupby(nfm.MOTIF_NAME_COL)):
    res = network_scan_build_single_tf(tf_mi_df)
    prior_matrix.append(res)
prior_matrix = pd.concat(prior_matrix, axis=1)
if extract_fasta and os.path.exists(extract_fasta):
    os.remove(extract_fasta)
prior_matrix = prior_matrix.astype(int)
prior_matrix.to_csv(os.path.join(path_out, 'prior.tsv'), sep='\t')
pd.DataFrame(prior_matrix.columns.values.reshape(-1, 1)).to_csv(os.path.join(path_out, 'tfs.tsv'), sep='\t', header=False, index=False)
