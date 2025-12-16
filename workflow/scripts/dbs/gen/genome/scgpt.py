from pathlib import Path
import sys


def ensure_checkpoint(model_dir):
    model_dir = Path(model_dir)
    vocab_file = model_dir / "vocab.json"
    args_file  = model_dir / "args.json"
    model_file = model_dir / "best_model.pt"
    if vocab_file.exists() and args_file.exists() and model_file.exists():
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
    gdown.download_folder(folder_url, output=str(model_dir), quiet=False, use_cookies=False)
    if not (vocab_file.exists() and args_file.exists() and model_file.exists()):
        raise FileNotFoundError(
            f"After attempted download, checkpoint files still missing in {model_dir}.\n"
            f"Please download manually from: {folder_url}"
        )


model_dir = sys.argv[1]
ensure_checkpoint(model_dir)
