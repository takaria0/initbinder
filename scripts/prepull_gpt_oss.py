#!/usr/bin/env python3
"""
Pre-pull (download) gpt-oss-120b weights to a specific directory, and optionally run a tiny test
generation via Transformers using the local folder (no network after download).

Usage examples:
  # Use a smaller model and custom dtype
  # Run with A100 sbatch job:
  sbatch -p gpu --gres=gpu:A100:1 --mem=64G -p ccl_lab_gpu --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=96G python prepull_gpt_oss.py --model-id openai/gpt-oss-120b --local-dir /pub/inagakit/.cache/huggingface --also-run
  
  python prepull_gpt_oss.py --model-id openai/gpt-oss-120b --local-dir /pub/inagakit/.cache/huggingface --also-run
"""

import os
import argparse
from pathlib import Path

def set_hf_paths(root: Path):
    """
    Override default Hugging Face paths to live under `root`.
    This affects both huggingface_hub and transformers caches.
    """
    root = '/pub/inagakit/.cache/huggingface' if root is None else root
    os.environ.setdefault("HF_HOME", str(root / ".hf_home"))
    os.environ.setdefault("HF_HUB_CACHE", str(root / ".hf_home" / "hub"))
    os.environ.setdefault("TRANSFORMERS_CACHE", str(root / ".transformers_cache"))
    # os.environ.setdefault("HF_HUB_DISABLE_TELEMETRY", "1")
    # Optional: speed up downloads if hf_transfer is installed
    # os.environ.setdefault("HF_HUB_ENABLE_HF_TRANSFER", "1")


def try_generate(model_id: Path, dtype: str, device_map: str, max_new_tokens: int):
    """
    Do a tiny test generation from the local folder using Transformers.
    """
    import torch
    from transformers import pipeline

    dtype_map = {
        "auto": "auto",
        "bfloat16": torch.bfloat16,
        "float16": torch.float16,
        "float32": torch.float32,
    }
    torch_dtype = dtype_map.get(dtype.lower(), "auto")

    print(f"[info] Initializing pipeline from {model_id} (dtype={dtype}, device_map={device_map}) …")
    pipe = pipeline(
        "text-generation",
        model=str(model_id),
        torch_dtype=torch_dtype,
        device_map=device_map,
    )

    # Minimal message to keep VRAM/compute tiny
    messages = [
        {"role": "user", "content": "Say 'hi'."},
    ]

    print("[info] Running a tiny generation to verify the model loads …")
    out = pipe(messages, max_new_tokens=max_new_tokens, do_sample=False, temperature=0.0)
    text = out[0].get("generated_text", out[0])
    if isinstance(text, list):
        # chat-style output
        last = text[-1]
        content = last["content"] if isinstance(last, dict) and "content" in last else str(last)
    else:
        content = str(text)
    print("[ok] Generation output (truncated):")
    print(content[:300])

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--model-id", default="openai/gpt-oss-120b", help="HF repo id")
    p.add_argument("--local-dir", required=True, help="Directory to store the model files")
    p.add_argument("--also-run", action="store_true", help="After download, run a tiny test generation")
    p.add_argument("--dtype", default="auto", choices=["auto", "bfloat16", "float16", "float32"], help="Torch dtype for test run")
    p.add_argument("--device-map", default="auto", help='Transformers device_map for test run (e.g., "auto", "cpu", "cuda")')
    p.add_argument("--max-new-tokens", type=int, default=1, help="Tokens for the tiny test run")
    args = p.parse_args()

    local_dir = Path(args.local_dir)

    # 1) Override default HF/Transformers caches under the same root (good hygiene)
    set_hf_paths(local_dir)

    # 2) Download the full snapshot directly into local_dir (no symlinks)
    # download_model(args.model_id, local_dir)
    print(f'[debug] Check if env is updated: HF_HOME={os.getenv("HF_HOME")}')

    # 3) Optional tiny generation to confirm it loads from disk
    if args.also_run:
        try:
            try_generate(args.model_id, args.dtype, args.device_map, args.max_new_tokens)
        except Exception as e:
            print(f"[warn] Test generation failed: {e}")
            print("       (Download is done; you can load later on a machine with enough memory.)")

if __name__ == "__main__":
    main()
