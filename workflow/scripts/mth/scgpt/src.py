"""
Attention-based GRN Inference with scGPT (single-pass attention capture, CPU ranking)

Usage (examples):
  python src.py \
      --data path/to/data.h5mu \
      --tfs TF.csv \
      --out output/grn.csv

  python src.py \
      --data path/to/data.h5mu \
      --tfs TF.csv \
      --out output/grn.csv \
      --attn_npy output/avg_attn.npy \
      --data_is_raw --n_hvg 2048 \
      --top_k_per_source 200 --min_score 0.8 \
      --batch_size 8 --device cuda:1
"""

import os
import warnings

os.environ["KMP_WARNINGS"] = "off"
warnings.filterwarnings(
    "ignore",
    message=r".*flash_attn is not installed.*",
    category=UserWarning
)

import json
import argparse
import logging
from pathlib import Path
from typing import List, Tuple, Dict
from tqdm import tqdm
import types
from contextlib import contextmanager

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F

# IO backends
import scanpy as sc
import mudata as mu
mu.set_options(pull_on_update=False)

from anndata import AnnData
from scipy.sparse import issparse

# scGPT imports
import scgpt as scg  # noqa: F401
from scgpt.preprocess import Preprocessor
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from scgpt.tokenizer import tokenize_and_pad_batch
from scgpt.model import TransformerModel
from scgpt.utils import set_seed

LOGGER = logging.getLogger("infer_grn_scgpt")

# ---------------------------
# Logging
# ---------------------------
def setup_logging(verbosity: int = 0, quiet: bool = False) -> None:
    if quiet:
        level = logging.WARNING
    else:
        level = logging.INFO if verbosity <= 0 else logging.DEBUG
    logging.basicConfig(
        level=level,
        format="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%H:%M:%S",
    )
    LOGGER.setLevel(level)

# ---------------------------
# Reproducibility
# ---------------------------
def set_seeds(seed=42):
    import random
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    set_seed(seed)

pad_token = "<pad>"
special_tokens = [pad_token, "<cls>", "<eoc>"]
pad_value = -2

# ---------------------------
# Parser
# ---------------------------
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Infer GRN using scGPT attention (single pass, CPU ranking).")
    p.add_argument("--data", required=True, help="Path to data (.h5ad or .h5mu).")
    p.add_argument("--model_dir", default=os.environ.get("SCGPT_MODEL_DIR", "./save/scGPT_human"),
                   help="Path to scGPT checkpoint dir (vocab.json, args.json, best_model.pt).")
    p.add_argument("--tfs", required=True, help="Path to headless one-column TF.csv (gene symbols).")
    p.add_argument("--out", required=True, help="Output GRN CSV (source,target,score).")
    p.add_argument("--attn_npy", default=None, help="Optional path to save the averaged attention (.npy).")

    # Preprocessing
    p.add_argument("--n_bins", type=int, default=51, help="Binning for expression.")
    p.add_argument("--n_hvg", type=int, default=2048, help="Number of HVGs (0 = no HVG filter).")
    p.add_argument("--data_is_raw", action="store_true",
                   help="Treat input .X as raw counts (apply log1p).")

    # Device / performance
    default_device = "cuda" if torch.cuda.is_available() else ("mps" if torch.backends.mps.is_available() else "cpu")
    p.add_argument("--device", default=default_device, help="Device: cuda | mps | cpu")
    p.add_argument("--batch_size", type=int, default=32, help="Batch size for encoder forward.")

    # GRN building
    p.add_argument("--top_k_per_source", type=int, default=1000, help="Top-K targets per TF.")
    p.add_argument("--min_score", type=float, default=0.75, help="Min score threshold for edges, 0=none, 1.0=max.")
    p.add_argument("--min_top_q", type=int, default=0, help="Ensure at least this many edges per TF.")
    p.add_argument("--drop_self_loops", action="store_true", help="Drop source->source edges.")
    p.add_argument("--restrict_targets_to_tfs", action="store_true",
                   help="Restrict candidate targets to TF list.")

    # UX
    p.add_argument("-v", "--verbose", action="count", default=0, help="Increase verbosity.")
    p.add_argument("-q", "--quiet", action="store_true", help="Only warnings/errors.")
    p.add_argument("--seed", type=int, default=42, help="Random seed.")
    return p

# ---------------------------
# Loaders
# ---------------------------
def load_adata(path: Path) -> AnnData:
    if path.suffix == ".h5mu":
        mdata = mu.read(path)
        return mdata.mod["rna"].copy()
    if path.suffix == ".h5ad":
        return sc.read(path)
    raise ValueError(f"Unsupported data format: {path.suffix}. Use .h5ad or .h5mu")

def ensure_checkpoint(model_dir: Path):
    vocab_file = model_dir / "vocab.json"
    args_file  = model_dir / "args.json"
    model_file = model_dir / "best_model.pt"
    if vocab_file.exists() and args_file.exists() and model_file.exists():
        LOGGER.debug("Found checkpoint files in %s", model_dir)
        return
    folder_url = "https://drive.google.com/drive/folders/1oWh_-ZRdhtoGQ2Fw24HP41FgLoomVo-y"
    model_dir.mkdir(parents=True, exist_ok=True)
    try:
        import gdown
    except ImportError:
        raise FileNotFoundError(
            f"Checkpoint not found in {model_dir} and gdown not installed.\n"
            f"Install with: pip install gdown\n"
            f"Or manually download from: {folder_url}"
        )
    LOGGER.info("Downloading pretrained checkpoint to %s ...", model_dir)
    gdown.download_folder(folder_url, output=str(model_dir), quiet=False, use_cookies=False)
    if not (vocab_file.exists() and args_file.exists() and model_file.exists()):
        raise FileNotFoundError(
            f"After attempted download, checkpoint files still missing in {model_dir}.\n"
            f"Please download manually from: {folder_url}"
        )

def load_checkpoint_and_vocab(model_dir_arg: str):
    model_dir = Path(model_dir_arg)
    ensure_checkpoint(model_dir)
    vocab_file = model_dir / "vocab.json"
    args_file  = model_dir / "args.json"
    model_file = model_dir / "best_model.pt"
    if not vocab_file.exists() or not args_file.exists() or not model_file.exists():
        raise FileNotFoundError(
            f"Could not find scGPT checkpoint files in {model_dir} "
            f"(expected vocab.json, args.json, best_model.pt)."
        )
    vocab = GeneVocab.from_file(vocab_file)
    with open(args_file, "r") as f:
        cfg = json.load(f)
    LOGGER.info("Loading model from %s (config: %s)", model_file, args_file)
    state = torch.load(model_file, map_location="cpu")
    return model_dir, vocab, cfg, state


# ---------------------------
# scGPT preprocessing
# ---------------------------
def _ensure_binned_layer(adata: AnnData, n_bins: int, data_is_raw: bool = False, subset_hvg: int = 0):
    pre = Preprocessor(
        use_key="X",
        filter_gene_by_counts=3,
        filter_cell_by_counts=True,
        normalize_total=1e4,
        result_normed_key="X_normed",
        log1p=data_is_raw,
        result_log1p_key="X_log1p",
        subset_hvg=subset_hvg if subset_hvg else False,
        hvg_flavor="seurat_v3" if data_is_raw else "cell_ranger",
        binning=n_bins,
        result_binned_key="X_binned",
    )
    pre(adata)
    return adata

def _intersect_with_vocab(adata: AnnData, vocab: GeneVocab) -> AnnData:
    adata.var["gene_name"] = adata.var_names
    adata.var["id_in_vocab"] = [1 if g in vocab else -1 for g in adata.var["gene_name"]]
    keep = (adata.var["id_in_vocab"] >= 0).values
    LOGGER.info("Intersected with vocab: kept %d / %d genes", int(keep.sum()), adata.n_vars)
    return adata[:, keep].copy()

def preprocess_and_intersect(
    adata: AnnData, vocab: GeneVocab, n_bins: int, data_is_raw: bool, n_hvg: int
) -> AnnData:
    adata = _ensure_binned_layer(adata, n_bins=n_bins, data_is_raw=data_is_raw, subset_hvg=n_hvg)
    return _intersect_with_vocab(adata, vocab)

def tokenize_adata(
    adata: AnnData, vocab: GeneVocab, pad_token: str, pad_value: int
) -> Tuple[np.ndarray, np.ndarray, torch.Tensor, List[str]]:
    Xb = adata.layers["X_binned"]
    all_counts = Xb.toarray() if issparse(Xb) else Xb
    gene_names = adata.var["gene_name"].astype(str).tolist()
    gene_ids = np.array(vocab(gene_names), dtype=int)
    tokenized_all = tokenize_and_pad_batch(
        all_counts,
        gene_ids,
        max_len=len(gene_names) + 1,
        vocab=vocab,
        pad_token=pad_token,
        pad_value=pad_value,
        append_cls=True,
        include_zero_gene=True,
    )
    all_gene_ids, all_values = tokenized_all["genes"], tokenized_all["values"]
    src_key_padding_mask = torch.as_tensor(all_gene_ids).eq(vocab[pad_token])
    return all_gene_ids, all_values, src_key_padding_mask, gene_names


# ---------------------------
# scGPT model
# ---------------------------
def build_model_from_cfg(vocab: GeneVocab, cfg: Dict, device: str):
    embsize   = cfg["embsize"]
    nhead     = cfg["nheads"]
    d_hid     = cfg["d_hid"]
    nlayers   = cfg["nlayers"]
    n_cls     = cfg.get("n_cls", 1)
    nlayers_cls = cfg.get("n_layers_cls", 1)
    pre_norm  = cfg.get("pre_norm", False)

    if pad_token not in vocab:
        vocab.append_token(pad_token)
    for s in special_tokens:
        if s not in vocab:
            vocab.append_token(s)

    model = TransformerModel(
        ntoken=len(vocab),
        d_model=embsize,
        nhead=nhead,
        d_hid=d_hid,
        nlayers=nlayers,
        nlayers_cls=nlayers_cls,
        n_cls=n_cls,
        vocab=vocab,
        dropout=0.0,
        pad_token=pad_token,
        pad_value=-2,
        domain_spec_batchnorm="batchnorm",
        use_fast_transformer=False,
        pre_norm=pre_norm,
    )

    loaded = False
    try:
        from scgpt.model import load_pretrained
        load_pretrained(model, state_dict=torch.load(
            Path(os.environ.get("SCGPT_MODEL_DIR", "./save/scGPT_human")) / "best_model.pt",
            map_location="cpu",
        ))
        loaded = True
    except Exception:
        pass
    if not loaded:
        sd = torch.load(Path(os.environ.get("SCGPT_MODEL_DIR", "./save/scGPT_human")) / "best_model.pt",
                        map_location="cpu")
        model.load_state_dict(sd, strict=False)

    model.eval().to(device)
    for p in model.parameters():
        p.requires_grad_(False)
    LOGGER.debug("Model built: %s", model)
    return model


# ---------------------------
# CUDA optimizations
# ---------------------------
def enable_fast_kernels_if_available():
    try:
        from torch.backends.cuda import sdp_kernel
        sdp_kernel(enable_flash=True, enable_mem_efficient=True, enable_math=False)
        torch.backends.cuda.matmul.allow_tf32 = True
        torch.backends.cudnn.allow_tf32 = True
        torch.set_float32_matmul_precision("high")
    except Exception as e:
        LOGGER.warning("Failed to enable fast SDPA: %s", e)

# ---------------------------
# Attention capture
# ---------------------------
@contextmanager
def capture_last_layer_attn(model, layer_idx: int = -1, average_attn_weights: bool = False):
    """
    Forces the last encoder layer's self-attention to return weights and captures them
    during a normal forward pass. Restores the module afterward.
    """
    layer = model.transformer_encoder.layers[layer_idx]
    mha = layer.self_attn
    captured = {}

    orig_forward = mha.forward

    def forward_patched(self, query, key, value, **kwargs):
        kwargs["need_weights"] = True
        kwargs["average_attn_weights"] = average_attn_weights
        out, w = orig_forward(query, key, value, **kwargs)
        captured["weights"] = w
        return out, w

    mha.forward = types.MethodType(forward_patched, mha)
    try:
        yield captured
    finally:
        mha.forward = orig_forward

def _to_BHSS(attn_w: torch.Tensor, B: int, S: int, num_heads: int) -> torch.Tensor:
    """
    Return shapes:
      - (B,H,S,S) unchanged
      - (B,S,S)   -> expand to (B,1,S,S) for a consistent API
      - (B*H,S,S) -> reshape to (B,H,S,S)
      - (S,B,H,S) -> permute to (B,H,S,S)
    """
    if attn_w.dim() == 4:
        if attn_w.shape[0] == B:
            return attn_w  # (B,H,S,S)
        return attn_w.permute(1, 2, 0, 3).contiguous()  # (S,B,H,S)->(B,H,S,S)
    if attn_w.dim() == 3:
        return attn_w.unsqueeze(1)  # (B,S,S)->(B,1,S,S)
    # fallback: (B*H,S,S)
    return attn_w.view(B, num_heads, S, S)

def aggregate_avg_attention(
    model: TransformerModel,
    all_gene_ids: np.ndarray,
    all_values: np.ndarray,
    src_key_padding_mask: torch.Tensor,
    batch_size: int,
    device: str,
) -> torch.Tensor:
    """
    Runs the single-pass encoder to capture last-layer attention, rank-normalizes
    rows and columns on CPU, and returns averaged (M x M) attention on CPU (float32).
    """
    layers = model.transformer_encoder.layers
    first = layers[0]
    batch_first = getattr(first, "batch_first", getattr(first.self_attn, "batch_first", False))
    norm_first  = getattr(first, "norm_first", False)

    N, S = all_gene_ids.shape[0], all_gene_ids.shape[1]
    M = S - 1

    # Accumulator for attention sum
    ACC_DEVICE = "cpu"
    sum_attn = torch.zeros((M, M), dtype=torch.float32, device=ACC_DEVICE)
    n_mats = 0

    LOGGER.info("TF32: matmul=%s, cudnn=%s, precision=%s",
                torch.backends.cuda.matmul.allow_tf32,
                torch.backends.cudnn.allow_tf32,
                torch.get_float32_matmul_precision(),)

    model.eval()
    with torch.inference_mode(), torch.cuda.amp.autocast(enabled=torch.cuda.is_available()):
        for i in tqdm(range(0, N, batch_size), desc="Extracting attention"):
            src  = torch.as_tensor(all_gene_ids[i:i+batch_size], dtype=torch.long,  device=device)
            vals = torch.as_tensor(all_values[i:i+batch_size],   dtype=torch.float, device=device)
            mask = src_key_padding_mask[i:i+batch_size].to(device)

            # encoder input
            total = model.encoder(src) + model.value_encoder(vals)  # (B,S,E)
            total = model.bn(total.permute(0, 2, 1)).permute(0, 2, 1)
            x = total if batch_first else total.transpose(0, 1)     # (B,S,E) or (S,B,E)

            with capture_last_layer_attn(model, layer_idx=len(layers)-1, average_attn_weights=True) as cap:
                # forward through all layers once
                for li, layer in enumerate(layers):
                    mha = layer.self_attn
                    dropout1 = getattr(layer, "dropout1", getattr(layer, "dropout", None))
                    dropout2 = getattr(layer, "dropout2", getattr(layer, "dropout", None))
                    act = layer.activation

                    sa_input = layer.norm1(x) if norm_first else x
                    attn_out, _ = mha(
                        sa_input, sa_input, sa_input,
                        key_padding_mask=mask,
                        need_weights=False  # last layer still captures due to patch
                    )

                    if li == len(layers) - 1:
                        # we only needed the weights; skip last FFN for a tiny win
                        x = sa_input
                    else:
                        if norm_first:
                            y = x + (dropout1(attn_out) if dropout1 is not None else attn_out)
                            ff_in  = layer.norm2(y)
                            ff_out = layer.linear2(layer.dropout(act(layer.linear1(ff_in))))
                            x      = y + (dropout2(ff_out) if dropout2 is not None else ff_out)
                        else:
                            y = layer.norm1(x + (dropout1(attn_out) if dropout1 is not None else attn_out))
                            ff_out = layer.linear2(layer.dropout(act(layer.linear1(y))))
                            x      = layer.norm2(y + (dropout2(ff_out) if dropout2 is not None else ff_out))

            # captured weights -> (B,H,S,S) -> mean over heads -> (B,S,S)
            mha_last = layers[-1].self_attn
            H = mha_last.num_heads
            w = cap["weights"]  # could be (B,S,S), (B,H,S,S) or (B*H,S,S)
            B = mask.shape[0]
            S_eff = src.shape[1] if batch_first else src.shape[0]
            w = _to_BHSS(w, B=B, S=S_eff, num_heads=H).mean(dim=1)  # (B,S,S)

            # Drop <CLS>, move to CPU
            A = w[:, 1:, 1:]
            A = A.to(torch.float16).to("cpu")
            cap.clear()
            del w
            if torch.cuda.is_available():
                torch.cuda.empty_cache()

            # ===== CPU rank-normalization =====
            A = A.to(dtype=torch.float32)

            # Row-wise ranks
            order = torch.argsort(A, dim=2)
            rank  = torch.argsort(order, dim=2).to(torch.float32) / M
            A = rank
            # Column-wise ranks
            order = torch.argsort(A.transpose(1, 2), dim=2)
            rank2  = torch.argsort(order, dim=2).to(torch.float32) / M
            A_rank = rank2.transpose(1, 2)  # back to (B,M,M)

            # Accumulate on CPU
            sum_attn.add_(A_rank.sum(dim=0))
            n_mats += A_rank.shape[0]
            del A, order, rank, rank2, A_rank

            if torch.cuda.is_available():
                torch.cuda.empty_cache()

    LOGGER.info("Aggregated attention across %d cell-batches", n_mats)
    avg_attn = (sum_attn / max(n_mats, 1))  # CPU float32 tensor
    return avg_attn


# ---------------------------
# GRN helpers
# ---------------------------
def _read_tf_names_case_insensitive(tfs_csv_path: str, attn_df: pd.DataFrame) -> List[str]:
    tf_list_upper = pd.read_csv(tfs_csv_path, header=None)[0].astype(str).str.upper().tolist()
    genes_upper_to_orig = {g.upper(): g for g in attn_df.columns.astype(str)}
    tf_present_upper = [u for u in tf_list_upper if u in genes_upper_to_orig]
    if not tf_present_upper:
        raise ValueError("None of the TFs in TF.csv are present in avg_attn_df columns.")
    return [genes_upper_to_orig[u] for u in tf_present_upper]

def _select_targets(attn_df: pd.DataFrame, tf_names: List[str], restrict_to_tfs: bool) -> pd.DataFrame:
    return attn_df.loc[:, tf_names] if restrict_to_tfs else attn_df

def _topk_threshold_indices(vals: np.ndarray, k_cap: int, min_score: float, min_top_q: int) -> np.ndarray:
    finite = np.isfinite(vals)
    if not finite.any():
        return np.array([], dtype=int)
    k = k_cap if (k_cap and k_cap > 0) else int(finite.sum())
    k = min(k, int(finite.sum()))
    if k == 0:
        return np.array([], dtype=int)
    top_idx = np.argpartition(-vals, kth=k-1)[:k]
    top_idx = top_idx[np.argsort(-vals[top_idx])]
    keep = top_idx[vals[top_idx] >= float(min_score)]
    if min_top_q and min_top_q > 0 and keep.size < min_top_q:
        fallback = top_idx[:min_top_q]
        seen = set(keep.tolist())
        merged = list(keep) + [i for i in fallback.tolist() if i not in seen]
        keep = np.array(merged, dtype=int)
    return keep

def _build_grn_from_attn(
    attn_df: pd.DataFrame,
    tf_names: List[str],
    top_k_per_source: int = 50,
    min_score: float = 0.01,
    min_top_q: int = 5,
    drop_self_loops: bool = True,
    restrict_targets_to_tfs: bool = False,
) -> pd.DataFrame:
    rows: List[Tuple[str, str, float]] = []
    attn_targets = _select_targets(attn_df, tf_names, restrict_targets_to_tfs)

    for src in tf_names:
        row = attn_targets.loc[src].astype(float).copy()
        if drop_self_loops and src in row.index:
            row.loc[src] = -np.inf
        row[~np.isfinite(row)] = -np.inf
        keep_idx = _topk_threshold_indices(
            vals=row.values,
            k_cap=top_k_per_source,
            min_score=min_score,
            min_top_q=min_top_q,
        )
        if keep_idx.size == 0:
            continue
        vals = row.values
        for j in keep_idx.tolist():
            rows.append((src, row.index[j], float(vals[j])))

    return pd.DataFrame(rows, columns=["source", "target", "score"])

def _write_grn_csv(grn: pd.DataFrame, out_path: str) -> None:
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    grn.to_csv(out_path, index=False)

def build_and_save_grn(
    avg_attn: np.ndarray,
    gene_names: List[str],
    tfs_path: str,
    out_csv: str,
    top_k: int,
    min_score: float,
    min_top_q: int,
    drop_self_loops: bool,
    restrict_targets_to_tfs: bool,
    attn_npy: str = None,
) -> pd.DataFrame:
    avg_attn_df = pd.DataFrame(avg_attn, columns=gene_names, index=gene_names)
    if attn_npy:
        Path(attn_npy).parent.mkdir(parents=True, exist_ok=True)
        np.save(attn_npy, avg_attn.astype(np.float32))
    tf_names = _read_tf_names_case_insensitive(tfs_path, avg_attn_df)
    grn = _build_grn_from_attn(
        avg_attn_df,
        tf_names,
        top_k_per_source=top_k,
        min_score=min_score,
        min_top_q=min_top_q,
        drop_self_loops=drop_self_loops,
        restrict_targets_to_tfs=restrict_targets_to_tfs,
    )
    _write_grn_csv(grn, out_csv)
    return grn

# ---------------------------
# Orchestration
# ---------------------------
def run(args: argparse.Namespace) -> None:
    setup_logging(args.verbose, args.quiet)
    set_seeds(args.seed)
    LOGGER.info("Device: %s", args.device)

    # Load checkpoint & vocab
    _model_dir, vocab, cfg, _state = load_checkpoint_and_vocab(args.model_dir)

    # Build model
    model = build_model_from_cfg(vocab, cfg, device=args.device)

    # Load & preprocess data
    data_path = Path(args.data)
    adata = load_adata(data_path)
    LOGGER.info("Loaded data: %s cells × %s genes", adata.n_obs, adata.n_vars)
    adata = preprocess_and_intersect(adata, vocab, args.n_bins, args.data_is_raw, args.n_hvg)

    # Tokenize
    all_gene_ids, all_values, src_key_padding_mask, gene_names = tokenize_adata(
        adata, vocab, pad_token, pad_value
    )

    # Performance toggles
    enable_fast_kernels_if_available()

    # Attention aggregation (single pass) with CPU ranking
    avg_attn_tensor = aggregate_avg_attention(
        model=model,
        all_gene_ids=all_gene_ids,
        all_values=all_values,
        src_key_padding_mask=src_key_padding_mask,
        batch_size=args.batch_size,
        device=args.device,
    )
    avg_attn = avg_attn_tensor.numpy()  # (M, M)

    # Build GRN (+ optional save of attention)
    grn = build_and_save_grn(
        avg_attn=avg_attn,
        gene_names=gene_names,
        tfs_path=args.tfs,
        out_csv=args.out,
        top_k=args.top_k_per_source,
        min_score=args.min_score,
        min_top_q=args.min_top_q,
        drop_self_loops=args.drop_self_loops,
        restrict_targets_to_tfs=args.restrict_targets_to_tfs,
        attn_npy=args.attn_npy,
    )

    LOGGER.info(f"Wrote GRN with {len(grn):,d} edges to {args.out}")
    LOGGER.info("Unique TF sources in GRN: %d", grn["source"].nunique() if not grn.empty else 0)

# ---------------------------
# CLI entry
# ---------------------------
def main():
    run(build_parser().parse_args())

if __name__ == "__main__":
    main()
