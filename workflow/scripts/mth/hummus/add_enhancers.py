import pandas as pd
import tqdm

def connect_genes_to_peaks_pandas_3(gene_peak_df: pd.DataFrame,
                                  peak_peak_df: pd.DataFrame,
                                  peak_cols=('Peak1','Peak2')) -> pd.DataFrame:
    """
    Connect each gene to peaks at distance 0, 1 or 2 through a peak–peak network
    using only pandas (no networkx).

    Parameters
    ----------
    gene_peak_df : pd.DataFrame
        Columns: ['gene','peak']
    peak_peak_df : pd.DataFrame
        Columns: [peak_cols[0], peak_cols[1], ...]
    peak_cols : tuple[str,str]
        Names of the two columns in the peak–peak network.

    Returns
    -------
    pd.DataFrame
        Columns: ['gene','peak','distance']
        Distance 0 = direct, 1 = one hop in peak network, 2 = two hops
    """

    peak_u, peak_v = peak_cols

    # 0-step (direct) connections
    direct_df = gene_peak_df[['gene','peak', 'distance']].drop_duplicates().copy()

    # 1-step: gene → peak (bipartite) → neighbor peaks (in peak network)
    step1_a = gene_peak_df.merge(peak_peak_df, left_on='peak', right_on=peak_u)[["gene", peak_v, 'distance']]
    step1_a = step1_a.rename(columns={peak_v: 'peak'})
    step1_b = gene_peak_df.merge(peak_peak_df, left_on='peak', right_on=peak_v)[["gene", peak_u, 'distance']]
    step1_b = step1_b.rename(columns={peak_u: 'peak'})
    step1_df = pd.concat([step1_a, step1_b], ignore_index=True).drop_duplicates()
    step1_df['distance_y'] = 1/10

    step1_df["distance"] = step1_df["distance"] * step1_df["distance_y"]
    step1_df = step1_df[['gene', 'peak', 'distance']]

    # 2-step: take step1 peaks and go one more hop
    step2_a = step1_df.merge(peak_peak_df, left_on='peak', right_on=peak_u)[['gene', peak_v]]
    step2_a = step2_a.rename(columns={peak_v: 'peak'})
    step2_b = step1_df.merge(peak_peak_df, left_on='peak', right_on=peak_v)[['gene', peak_u]]
    step2_b = step2_b.rename(columns={peak_u: 'peak'})
    step2_df = pd.concat([step2_a, step2_b], ignore_index=True).drop_duplicates()
    step2_df['distance'] = 1/100

    # 3-step: take step2 peaks and go one more hop
    step3_a = step2_df.merge(peak_peak_df, left_on='peak',
                            right_on=peak_u)[['gene', peak_v]]
    step3_a = step3_a.rename(columns={peak_v: 'peak'})
    step3_b = step2_df.merge(peak_peak_df, left_on='peak',
                            right_on=peak_v)[['gene', peak_u]]
    step3_b = step3_b.rename(columns={peak_u: 'peak'})
    step3_df = pd.concat([step3_a, step3_b], ignore_index=True).drop_duplicates()
    step3_df['distance'] = 1/1000

    # Combine all, keep shortest distance for duplicates
    all_steps = pd.concat([direct_df, step1_df, step2_df, step3_df], ignore_index=True)

    # Combine all, keep shortest distance for duplicates
    result = (all_steps.sort_values(['gene','peak','distance'])
                      .drop_duplicates(subset=['gene','peak'], keep='first')
                      .reset_index(drop=True))

    return result


def connect_genes_to_peaks_via_genes_4(gene_peak_df: pd.DataFrame,
                                       gene_gene_df: pd.DataFrame,
                                       gene_cols=('Gene1','Gene2')) -> pd.DataFrame:
    """
    Connect each source gene to peaks at distance 0..4 where distance counts hops
    through the gene–gene network (no networkx, pure pandas).

    Parameters
    ----------
    gene_peak_df : pd.DataFrame
        Columns: ['gene','peak']  (bipartite edges)
    gene_gene_df : pd.DataFrame
        Columns: [gene_cols[0], gene_cols[1], ...]  (gene–gene edges)
    gene_cols : tuple[str,str]
        Names of the two columns in the gene–gene network.

    Returns
    -------
    pd.DataFrame
        Columns: ['gene','peak','distance']
        Distance 0 = direct GP edge,
        Distance k>=1 = k hops in GG from the source gene, then attach peaks of reached genes.
    """
    g_u, g_v = gene_cols

    # Normalize
    gp = gene_peak_df[['gene','peak']].dropna().astype(str).drop_duplicates().copy()
    gg = gene_gene_df[[g_u, g_v]].dropna().astype(str).drop_duplicates().copy()

    # Seed set of source genes (use the genes present in gp; extend with gg if you like)
    sources = pd.DataFrame({'gene': gp['gene'].unique()})

    # Build undirected GG neighbor table
    gg_a = gg.rename(columns={g_u: 'src', g_v: 'dst'})[['src','dst']]
    gg_b = gg.rename(columns={g_u: 'dst', g_v: 'src'})[['src','dst']]
    neighbors = pd.concat([gg_a, gg_b], ignore_index=True).drop_duplicates()

    # ---- distance 0: direct gene→peak
    direct_df = gp.copy()
    direct_df['distance'] = 1

    # Helper: attach peaks of the reached genes to produce (source gene, peak)
    def attach_peaks(reached_gene_pairs: pd.DataFrame) -> pd.DataFrame:
        # reached_gene_pairs columns: ['gene','reached_gene']
        att = (reached_gene_pairs.merge(gp, left_on='reached_gene', right_on='gene', how='left')
                               .dropna(subset=['peak']))
        att = att.rename(columns={'gene_x':'gene'}).loc[:, ['gene','peak']].drop_duplicates()
        return att

    # ---- step1: 1 GG hop
    step1_genes_a = sources.merge(gg, left_on='gene', right_on=g_u)[['gene', g_v]].rename(columns={g_v:'reached_gene'})
    step1_genes_b = sources.merge(gg, left_on='gene', right_on=g_v)[['gene', g_u]].rename(columns={g_u:'reached_gene'})
    step1_genes = pd.concat([step1_genes_a, step1_genes_b], ignore_index=True).drop_duplicates()
    step1_df = attach_peaks(step1_genes)
    step1_df['distance'] = 1/100

    # Combine all, keep shortest distance per (gene, peak)
    all_steps = pd.concat([direct_df, step1_df], ignore_index=True)
    result = (all_steps.sort_values(['gene','peak','distance'])
                      .drop_duplicates(subset=['gene','peak'], keep='first')
                      .reset_index(drop=True))
    return result


def topk_peaks_per_tf(
    tf_peak_df: pd.DataFrame,         # columns: ['tf','peak','score']
    gene_peak_df: pd.DataFrame,       # columns: ['gene','peak','score']
    k: int = 3,                       # keep k peaks per (tf, gene)
    top_m_per_tf: int | None = None,  # optional: prefilter top-m peaks per TF before merging
) -> pd.DataFrame:
    """
    For each TF, keep the top-k peaks per gene by product score = score_tf_peak * score_gene_peak.
    Processes TFs one-by-one to avoid a huge global merge.

    Returns columns: ['tf','peak','gene','score']
    """
    # keep only needed columns
    tf_small = tf_peak_df[['tf','peak','score']].dropna().copy()
    gp_small = gene_peak_df[['gene','peak','score']].dropna().copy()

    # ensure consistent dtypes for 'peak'
    if tf_small['peak'].dtype != gp_small['peak'].dtype:
        tf_small['peak'] = tf_small['peak'].astype(str)
        gp_small['peak'] = gp_small['peak'].astype(str)

    # optional: limit to shared peaks to shrink work
    common_peaks = pd.Index(tf_small['peak']).intersection(gp_small['peak'])
    tf_small = tf_small[tf_small['peak'].isin(common_peaks)]
    gp_small = gp_small[gp_small['peak'].isin(common_peaks)]

    # pre-sort gp_small by (peak, score desc) so groupby operations are fast
    gp_small = gp_small.sort_values(['peak','score'], ascending=[True, False])

    out_parts = []

    # iterate per TF (keeps memory bounded)
    for tf, block in tqdm.tqdm(tf_small.groupby('tf', sort=False)):
        tf_block = block

        # optional prefilter: take top-m strongest peaks for this TF
        if top_m_per_tf is not None and len(tf_block) > top_m_per_tf:
            tf_block = tf_block.nlargest(top_m_per_tf, 'score')

        # merge this TF's peaks with genes on the same peaks
        m = tf_block.merge(
            gp_small, on='peak', how='inner', suffixes=('_tf','_gene')
        )
        if m.empty:
            continue

        # per-peak contribution = product
        m['score'] = m['score_tf'] * m['score_gene']
        m = m[['tf','peak','gene','score']]

        # keep top-k peaks per gene for this TF
        m = (m.sort_values(['gene','score'], ascending=[True, False])
               .groupby('gene', as_index=False, group_keys=False)
               .head(k))

        out_parts.append(m)

    if not out_parts:
        return pd.DataFrame(columns=['tf','peak','gene','score'])

    out = pd.concat(out_parts, ignore_index=True)

    # optional final tidy: if a (tf,gene,peak) appears multiple times, dedupe on best score
    out = (out.sort_values(['tf','gene','peak','score'], ascending=[True, True, True, False])
             .drop_duplicates(subset=['tf','gene','peak'], keep='first')
             .reset_index(drop=True))

    return out

import os
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-b', '--binding_regions', required=True)
parser.add_argument('-t', '--atac_rna', required=True)
parser.add_argument('-r', '--rna_layer', required=True)
parser.add_argument('-a', '--atac_layer', required=True)
parser.add_argument('-i', '--intermediar_grn', required=True)
parser.add_argument('-o', '--final_grn', required=True)
parser.add_argument('-g', '--genome_annotations', required=True)
args = vars(parser.parse_args())

# Parameters
num_cres = 10  # number of enhancers to keep per (TF, gene

# File paths
tf_atac_path = args['binding_regions']
atac_rna_path = args['atac_rna']
rna_layer_path = args['rna_layer']
atac_layer_path = args['atac_layer']
intermediar_grn_path = args['intermediar_grn']
final_grn_path = args['final_grn']
genome_annotations_path = args['genome_annotations']


# Binding regions
tf_atac = pd.read_csv(tf_atac_path)[["tf", "peak", "score"]]
tf_atac.tf = tf_atac.tf.astype('category')
tf_atac.peak = tf_atac.peak.astype('category')

# Direct promoter links
atac_rna = pd.read_csv(atac_rna_path, header=None, sep='\t')
atac_rna.columns = ['gene', 'peak']

# Intermediate layers to connect enhancers to genes
atac_layer = pd.read_csv(atac_layer_path)
rna_layer = pd.read_csv(rna_layer_path)
rna_layer.columns = ['gene1', 'gene2', 'score']
rna_layer = rna_layer.iloc[:50000, :]

# Connect peaks regulating the same gene
atac_layer_from_rna = atac_rna.rename(columns={'peak':'Peak1'}).merge(atac_rna.rename(columns={'peak':'Peak2'}), on='gene')
# Add them to the atac layer
atac_layer = pd.concat([
    atac_layer[['Peak1', 'Peak2']],
    atac_layer_from_rna[['Peak1', 'Peak2']]
], ignore_index=True).drop_duplicates()

# Connect peaks to genes via the intermediate layers

## Connect genes to peaks via genes
gene_peak_df = connect_genes_to_peaks_via_genes_4(
    atac_rna, rna_layer, gene_cols=["gene1", "gene2"])
gene_peak_df.gene = gene_peak_df.gene.astype('category')
gene_peak_df.peak = gene_peak_df.peak.astype('category')

## Connect genes to peaks via peaks
gene_peak_df = connect_genes_to_peaks_pandas_3(
    gene_peak_df, atac_layer, peak_cols=["Peak1", "Peak2"])
gene_peak_df.columns = ['gene', 'peak', 'score']

# filter to keep peaks only on the gene chromosome
gene_mapping = pd.read_csv(
    os.path.join(genome_annotations_path)
)
gene_mapping = gene_mapping[['gene_name', 'seqnames']].drop_duplicates()
gene_mapping.set_index('gene_name', inplace=True)
print(f"Number of genes with chromosome info: {gene_mapping.head()}")

gene_mapping = gene_mapping['seqnames'].to_dict()
gene_peak_df = gene_peak_df[gene_peak_df.peak.str.split('-').str[0] == gene_peak_df.gene.map(gene_mapping)]

del atac_layer, atac_layer_from_rna, rna_layer, atac_rna

common_peaks = tf_atac.peak.unique()[pd.Series(tf_atac.peak.unique()).isin(gene_peak_df.peak.unique()).values]

# Merge TF-peak and gene-peak on peaks, keep top-k peaks per (TF, gene)
topk_peaks_per_tf_df = topk_peaks_per_tf(
    tf_atac[tf_atac.peak.isin(common_peaks)],
    gene_peak_df[gene_peak_df.peak.isin(common_peaks)],
    k=num_cres,)
topk_peaks_per_tf_df.columns = ["TF", "CRE", "gene", "enhancer_score"]

# Load TF - gene scores
grn = pd.read_csv(intermediar_grn_path)
grn.rename(columns = {"source": "TF", "target": "gene", "score": "score"}, inplace = True)

# Filter to TFs and genes present in the triplets
grn = grn[grn.TF.isin(topk_peaks_per_tf_df.TF.unique()) & grn.gene.isin(topk_peaks_per_tf_df.gene.unique())]

# Merge with GRN to keep only high-confidence TF-gene pairs
final = topk_peaks_per_tf_df.merge(grn, on=["TF", "gene"], how="inner")
final = final[["TF", "CRE", "gene", "score"]]

final.to_csv(
    final_grn_path,
    sep="\t", index=False, header=True)
