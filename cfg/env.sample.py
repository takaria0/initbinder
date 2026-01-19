import os
from pathlib import Path

# Copy this file to cfg/env.py and fill in real values.
# Keep secrets out of git.

# --- PyMOL bundle destination (overridable via env) ---
_default_pymol_dest = Path.cwd() / "pymol_bundles"
RFA_LOCAL_PYMOL_DEST = str(Path(os.getenv("RFA_LOCAL_PYMOL_DEST", str(_default_pymol_dest))).expanduser())
RFA_LOCAL_PYMOL_SSH_OPTS = os.getenv(
    "RFA_LOCAL_PYMOL_SSH_OPTS",
    "-o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null",
)

# --- LLM configuration knobs ---
USE_LLM = os.getenv("USE_LLM", "false").strip().lower() == "true"
LLM_PROVIDER = os.getenv("LLM_PROVIDER", "openai").strip().lower()
MODEL = os.getenv("MODEL", "gpt-4.1-2025-04-14")

OPENAI_API_KEY = os.getenv("OPENAI_API_KEY", "YOUR_OPENAI_KEY")
GROQ_API_KEY = os.getenv("GROQ_API_KEY", "YOUR_GROQ_KEY")
GOOGLE_API_KEY = os.getenv("GOOGLE_API_KEY", "YOUR_GOOGLE_KEY")

# Optional knobs for local Transformers backends
TORCH_DTYPE = os.getenv("TORCH_DTYPE", "auto")
DEVICE_MAP = os.getenv("DEVICE_MAP", "auto")
MAX_NEW_TOKENS = int(os.getenv("MAX_NEW_TOKENS", "8096"))
REASONING_EFFORT = os.getenv("REASONING_EFFORT", "high")

# --- PyMOL mode ---
RFA_PYMOL_MODE = os.getenv("RFA_PYMOL_MODE", "bundle")
RFA_PYMOL_REMOTE_HOST = os.getenv("RFA_PYMOL_REMOTE_HOST", "localhost")
RFA_PYMOL_REMOTE_PORT = int(os.getenv("RFA_PYMOL_REMOTE_PORT", "9124"))

# --- Hugging Face cache roots ---
HF_ROOT = Path(os.getenv("HF_ROOT", str(Path.home() / ".cache" / "huggingface")))
os.environ.setdefault("HF_HOME", str(HF_ROOT / ".hf_home"))
os.environ.setdefault("HF_HUB_CACHE", str(HF_ROOT / ".hf_home" / "hub"))
os.environ.setdefault("TRANSFORMERS_CACHE", str(HF_ROOT / ".transformers_cache"))
