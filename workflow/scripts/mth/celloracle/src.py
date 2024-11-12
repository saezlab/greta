import matplotlib.pyplot as plt  # Import else compiling error
import pandas as pd
import numpy as np
from celloracle import motif_analysis as ma
import celloracle as co
import mudata as mu
import os
import argparse
from genomepy import Genome, install_genome, config
from gimmemotifs.motif import default_motifs


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-a','--path_data', required=True)
parser.add_argument('-b','--all_peaks', required=True)
parser.add_argument('-c','--connections', required=True)
parser.add_argument('-d','--organism', required=True)
parser.add_argument('-e','--thr', required=True)
parser.add_argument('-f','--fpr', required=True)
parser.add_argument('-g','--blen', required=True)
parser.add_argument('-i','--tfb_thr', required=True)
parser.add_argument('-j','--alpha', required=True)
parser.add_argument('-k','--pthr', required=True)
parser.add_argument('-l','--top_n', required=True)
parser.add_argument('-m','--knn', required=True)
parser.add_argument('-n','--path_out', required=True)
args = vars(parser.parse_args())

path_data = args['path_data']
path_all_peaks = args['all_peaks']
path_connections = args['connections']
org = args['organism']
thr_coaccess = float(args['thr'])
fpr = float(args['fpr'])
blen = int(args['blen'])
tfb_thr = float(args['tfb_thr'])
alpha = float(args['alpha'])
pthr = float(args['pthr'])
top_n = int(args['top_n'])
k = int(args['knn'])
path_out = args['path_out']

n_jobs = 32

if __name__ == '__main__':
    # Load scATAC-seq peak list
    peaks = pd.read_csv(path_all_peaks, index_col=0).x.values.astype('U')
    peaks = np.char.replace(peaks, '-', '_')
    
    # Load Cicero coaccessibility scores
    cicero_connections = pd.read_csv(path_connections, index_col=0)
    cicero_connections['Peak1'] = np.char.replace(cicero_connections['Peak1'].values.astype('U'), '-', '_')
    cicero_connections['Peak2'] = np.char.replace(cicero_connections['Peak2'].values.astype('U'), '-', '_')
    
    # Extract tss information
    tss_annotated = ma.get_tss_info(
        peak_str_list=peaks,
        ref_genome=org
    )
    
    # Integrate
    p2g = ma.integrate_tss_peak_with_cicero(
        tss_peak=tss_annotated,
        cicero_connections=cicero_connections
    )
    
    # Process
    p2g = p2g[p2g['coaccess'] >= thr_coaccess]
    p2g['peak_id'] = p2g['peak_id'].str.replace('_', '-')
    p2g = p2g.rename(columns={'peak_id': 'cre', 'gene_short_name': 'gene', 'coaccess': 'score'})
    p2g = p2g.sort_values(['cre', 'score'], ascending=[True, False])
    
    #################################################
    
    # Load annotated peak data.
    p2g['cre'] = p2g['cre'].str.replace('-', '_')
    
    def decompose_chrstr(peak_str):
        """
        Args:
            peak_str (str): peak_str. e.g. 'chr1_3094484_3095479'
    
        Returns:
            tuple: chromosome name, start position, end position
        """
    
        *chr_, start, end = peak_str.split("_")
        chr_ = "_".join(chr_)
        return chr_, start, end
    
    def check_peak_format(peaks_df, gname, gdir):
        """
        Check peak format.
         (1) Check chromosome name.
         (2) Check peak size (length) and remove sort DNA sequences (<5bp)
    
        """
    
        df = peaks_df.copy()
        df = df.rename(columns={'cre': 'peak_id', 'gene': 'gene_short_name'})
        n_peaks_before = df.shape[0]
    
        # Decompose peaks and make df
        decomposed = [decompose_chrstr(peak_str) for peak_str in df["peak_id"]]
        df_decomposed = pd.DataFrame(np.array(decomposed), index=peaks_df.index)
        df_decomposed.columns = ["chr", "start", "end"]
        df_decomposed["start"] = df_decomposed["start"].astype(int)
        df_decomposed["end"] = df_decomposed["end"].astype(int)
    
        # Load genome data
        genome_data = Genome(gname, genomes_dir=gdir)
        all_chr_list = list(genome_data.keys())
    
        # DNA length check
        lengths = np.abs(df_decomposed["end"] - df_decomposed["start"])
    
        # Filter peaks with invalid chromosome name
        n_threshold = 5
        df = df[(lengths >= n_threshold) & df_decomposed.chr.isin(all_chr_list)]
    
        # DNA length check
        lengths = np.abs(df_decomposed["end"] - df_decomposed["start"])
    
        # Data counting
        n_invalid_length = len(lengths[lengths < n_threshold])
        n_peaks_invalid_chr = n_peaks_before - df_decomposed.chr.isin(all_chr_list).sum()
        n_peaks_after = df.shape[0]
    
        # Print
        print("Peaks before filtering: ", n_peaks_before)
        print("Peaks with invalid chr_name: ", n_peaks_invalid_chr)
        print("Peaks with invalid length: ", n_invalid_length)
        print("Peaks after filtering: ", n_peaks_after)
    
        return df
    
    # Format and delete peaks
    gdir = f'dbs/{org}/gen/genome/celloracle/'
    p2g = check_peak_format(p2g, gname=org, gdir=gdir)
    
    # Instantiate TFinfo object
    tfi = ma.TFinfo(
        peak_data_frame=p2g,
        ref_genome=org,
        genomes_dir=gdir,
    )
    
    # Update config
    config.config.config['genomes_dir'] = gdir
    
    # Scan
    tfi.scan(
        background_length=blen,
        fpr=fpr,
        motifs=None,
        verbose=True,
        n_cpus=n_jobs,
    )
    
    # Do filtering
    tfi.filter_motifs_by_score(threshold=tfb_thr)
    
    # Extract filtered TF predictions
    tfb = tfi.scanned_filtered[["seqname", "motif_id", "score"]].copy()
    tfb['motif_values'] = tfb['motif_id'].map(tfi.dic_motif2TFs)
    tfb = tfb.explode('motif_values')
    tfb = tfb.groupby(['seqname', 'motif_values'])['score'].max().reset_index()
    tfb = tfb[['seqname', 'motif_values', 'score']].dropna()
    tfb = tfb.reset_index(drop=True).rename(columns={'seqname': 'cre', 'motif_values': 'tf'})
    tfb['cre'] = tfb['cre'].str.replace('_', '-')
    tfb = tfb.sort_values(['cre', 'score'], ascending=[True, False])
    
    ########################################################
    
    # Process base GRN
    tfb['score'] = 1
    p2g = p2g.rename(columns={'peak_id': 'cre', 'gene_short_name': 'gene'})
    p2g['cre'] = p2g['cre'].str.replace('_', '-')
    p2g = p2g[['cre', 'gene']]
    base_grn = pd.merge(
        p2g,
        tfb
        .pivot(index='cre', columns='tf')
        .fillna(0)
        .droplevel(0, axis=1)
        .reset_index()
    )
    base_grn = base_grn.rename(columns={'cre': 'peak_id', 'gene': 'gene_short_name'})
    base_grn['peak_id'] = base_grn['peak_id'].str.replace('-', '_')
    
    # Init oracle object
    # Extract raw counts data and assign labels
    mdata = mu.read(path_data)
    adata = mdata.mod['rna'].copy()
    adata.layers['lognorm'] = adata.X.copy()
    adata.X = adata.layers['counts'].copy()
    adata.obs['cluster'] = 'cluster'
    adata.obsm['X_pca'] = mdata.obsm['X_spectral']
    
    # Instantiate Oracle object
    oracle = co.Oracle()
    oracle.import_anndata_as_raw_count(
        adata=adata,
        cluster_column_name="cluster",
        embedding_name="X_pca"
    )
    
    # Compute PCA and select top pcs
    oracle.perform_PCA()
    n_comps = np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_))>0.002))[0][0]
    n_comps = min(n_comps, 50)
    
    # Run imputation
    oracle.knn_imputation(
        n_pca_dims=n_comps,
        k=k,
        balanced=True,
        b_sight=k*8,
        b_maxl=k*4,
        n_jobs=n_jobs,
    )
    oracle.import_TF_data(TF_info_matrix=base_grn)
    
    # Model TF ~ G
    print('Modeling GRN...')
    links = oracle.get_links(
        cluster_name_for_GRN_unit="cluster",
        alpha=alpha,
        n_jobs=n_jobs,
    )
    print('Modeling Done!')
    print('Filtering links...')
    links.filter_links(
        p=pthr,
        weight="coef_abs",
        threshold_number=top_n
    )
    print('Filtering done!')
    
    # Extract grn
    grn = links.filtered_links['cluster'].dropna()[['source', 'target', 'coef_mean', 'p']]
    grn = grn.rename(columns={'coef_mean': 'score', 'p': 'pval'})
    grn = grn.sort_values(['source', 'target', 'pval'])
    
    # Write
    grn.to_csv(path_out, index=False)
    
    print('Done')
    os._exit(0)  # Add this else it gets stuck