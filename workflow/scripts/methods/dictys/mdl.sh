#!/bin/bash


# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --output_d) output_d="$2"; shift ;;
        --pre_path) pre_path="$2"; shift ;;
        --p2g_path) p2g_path="$2"; shift ;;
        --tfb_path) tfb_path="$2"; shift ;;
        --annot) annot="$2"; shift ;;
        --distance) distance="$2"; shift ;;
        --n_p2g_links) n_p2g_links="$2"; shift ;;
        --threads) threads="$2"; shift ;;
        --device) device="$2"; shift ;;
        --out_path) out_path="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done && \
if [ $(wc -l < $p2g_path) -eq 1 ] || [ $(wc -l < $tfb_path) -eq 1 ]; then
    awk 'BEGIN {{ print "source,target,score,pval" }}' > $output_out
    exit 0
fi && \
mkdir -p "$output_d" && \
python -c "import torch; print('Cuda enabled:', torch.cuda.is_available())" && \
python -c "import pandas as pd, numpy as np, mudata as mu, sys, os; \
tfb = pd.read_csv(sys.argv[1]); \
tfb['cre'] = tfb['cre'].str.replace('-', ':'); \
peaks = tfb['cre'].unique(); \
rna = mu.read(os.path.join(sys.argv[2], 'mod', 'rna')); \
pd.DataFrame(np.zeros((peaks.size, 1)), index=peaks, columns=['placeholder']).to_csv(sys.argv[3], sep='\t', compression='gzip'); \
rna.to_df().T.to_csv(sys.argv[4], sep='\t', compression='gzip'); \
tfb.rename(columns={'cre': 'loc', 'tf': 'TF'})[['TF', 'loc', 'score']].to_csv(sys.argv[5], sep='\t', index=False)" \
$tfb_path $pre_path $output_d/peaks.tsv.gz $output_d/expr.tsv.gz $output_d/tfb.tsv.gz && \
python -m dictys chromatin tssdist --cut $distance $output_d/expr.tsv.gz $output_d/peaks.tsv.gz $annot $output_d/tssdist.tsv.gz && \
echo 'Finished tssdist' && \
python -m dictys chromatin linking $output_d/tfb.tsv.gz $output_d/tssdist.tsv.gz $output_d/linking.tsv.gz && \
echo 'Finished chromatin linking' && \
python -m dictys chromatin binlinking $output_d/linking.tsv.gz $output_d/binlinking.tsv.gz $n_p2g_links && \
echo 'Finished chromatin binlinking' && \
python -m dictys network reconstruct --device $device --nth $threads $output_d/expr.tsv.gz $output_d/binlinking.tsv.gz $output_d/net_weight.tsv.gz $output_d/net_meanvar.tsv.gz $output_d/net_covfactor.tsv.gz $output_d/net_loss.tsv.gz $output_d/net_stats.tsv.gz && \
echo 'Finished network reconstruct' && \
python -m dictys network normalize --nth $threads $output_d/net_weight.tsv.gz $output_d/net_meanvar.tsv.gz $output_d/net_covfactor.tsv.gz $output_d/net_nweight.tsv.gz && \
echo 'Finished network normalize' && \
python -c "import pandas as pd, numpy as np, sys, os; \
weights = pd.read_csv(sys.argv[1], sep='\t', index_col=0); \
mask = pd.read_csv(sys.argv[2], sep='\t', index_col=0); \
mask = mask.loc[weights.index, weights.columns]; \
df = [(weights.index[i], weights.columns[j], weights.iloc[i, j]) for i in np.arange(weights.shape[0]) for j in np.arange(weights.shape[1]) if mask.iloc[i, j]]; \
df = np.array(df); \
df = pd.DataFrame(df, columns=['source', 'target', 'score']); \
df['pval'] = 0.01; \
df.to_csv(sys.argv[3], index=False)" $output_d/net_nweight.tsv.gz $output_d/binlinking.tsv.gz $out_path
