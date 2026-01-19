InitBinder Web UI (Bulk)

Quickstart (local)
1) Create a virtualenv and install deps:
```bash
python3 -m venv .venv
. .venv/bin/activate
pip install -r requirements-webapp.txt
```

2) (Optional) Copy the default config and edit paths if needed:
```bash
cp cfg/webapp.yaml cfg/webapp.local.yaml
```
Use `cfg/webapp.local.yaml` if you want to keep local overrides. It is loaded after
`cfg/webapp.yaml` and is the best place to put machine-specific paths.

Config precedence
- `cfg/webapp.yaml` (defaults)
- `cfg/webapp.local.yaml` (local overrides)
- `INITBINDER_*` environment variables (highest priority)

Common overrides
- `INITBINDER_PROJECT_ROOT`
- `INITBINDER_TARGETS_DIR`
- `INITBINDER_REMOTE_ROOT`
- `INITBINDER_TARGET_ROOT`
- `INITBINDER_SSH_CONTROL_PATH`

3) Run the web server:
```bash
uvicorn webapp.main:app --reload --port 8000
```

4) Open the Bulk UI:
- http://localhost:8000/bulk

Notes
- The Bulk UI works without PyMOL or cluster access; those features are optional.
- Use any TSV/CSV pasted into the Bulk page; a small sample can be taken from `targets_catalog/*.tsv`.
- If you need to change paths/log locations, edit `cfg/webapp.local.yaml` or set env vars.
- If you use LLM features, copy `cfg/env.sample.py` to `cfg/env.py` and fill in keys (keep secrets out of git).
