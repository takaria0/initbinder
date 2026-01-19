import os

# Copy this file to cfg/env.py and fill in real values.
# Keep secrets out of git.

# --- LLM configuration knobs ---
USE_LLM = os.getenv("USE_LLM", "false").strip().lower() == "true"
LLM_PROVIDER = os.getenv("LLM_PROVIDER", "openai").strip().lower()
MODEL = os.getenv("MODEL", "gpt-4.1-2025-04-14")

OPENAI_API_KEY = os.getenv("OPENAI_API_KEY", "YOUR_OPENAI_KEY")
GROQ_API_KEY = os.getenv("GROQ_API_KEY", "YOUR_GROQ_KEY")
GOOGLE_API_KEY = os.getenv("GOOGLE_API_KEY", "YOUR_GOOGLE_KEY")

# Optional LLM limits
MAX_NEW_TOKENS = int(os.getenv("MAX_NEW_TOKENS", "8096"))
