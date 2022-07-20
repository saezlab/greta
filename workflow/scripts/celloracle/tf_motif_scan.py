import pandas as pd
import numpy as np
import celloracle as co
from celloracle import motif_analysis as ma
from celloracle.utility import save_as_pickled_object
from genomepy import Genome
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-p','--peaks', required=True)
parser.add_argument('-o','--organism', required=True)
parser.add_argument('-f','--fpr', required=True)
parser.add_argument('-t','--tfinfo', required=True)
args = vars(parser.parse_args())

path_peaks = args['peaks']
organism = args['organism']
fpr = float(args['fpr'])
path_tfinfo = args['tfinfo']

# Determine genome
if organism == 'human':
    genome = 'hg38'
elif organism == 'mouse':
    genome = 'mm10'

# Load annotated peak data.
peaks = pd.read_csv(path_peaks, index_col=0)

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


def check_peak_format(peaks_df, ref_genome):
    """
    Check peak format.
     (1) Check chromosome name.
     (2) Check peak size (length) and remove sort DNA sequences (<5bp)

    """

    df = peaks_df.copy()

    n_peaks_before = df.shape[0]

    # Decompose peaks and make df
    decomposed = [decompose_chrstr(peak_str) for peak_str in df["peak_id"]]
    df_decomposed = pd.DataFrame(np.array(decomposed), index=peaks_df.index)
    df_decomposed.columns = ["chr", "start", "end"]
    df_decomposed["start"] = df_decomposed["start"].astype(int)
    df_decomposed["end"] = df_decomposed["end"].astype(int)

    # Load genome data
    genome_data = Genome(ref_genome)
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


    #
    print("Peaks before filtering: ", n_peaks_before)
    print("Peaks with invalid chr_name: ", n_peaks_invalid_chr)
    print("Peaks with invalid length: ", n_invalid_length)
    print("Peaks after filtering: ", n_peaks_after)

    return df

# Check if right genome
peaks = check_peak_format(peaks, genome)

# Instantiate TFinfo object
tfi = ma.TFinfo(peak_data_frame=peaks,
                ref_genome=genome)

# Scan motifs. !!CAUTION!! This step may take several hours if you have many peaks!
tfi.scan(fpr=fpr,
         motifs=None,  # If you enter None, default motifs will be loaded.
         verbose=True)

# Save
tfi.to_hdf5(file_path=path_tfinfo)

