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
Use `cfg/webapp.local.yaml` if you want to keep local overrides.

3) Run the web server:
```bash
uvicorn webapp.main:app --reload --port 8000
```

4) Open the Bulk UI:
- http://localhost:8000/bulk

Notes
- The Bulk UI works without PyMOL or cluster access; those features are optional.
- Use any TSV/CSV pasted into the Bulk page; a small sample can be taken from `targets_catalog/*.tsv`.
- If you need to change paths/log locations, edit `cfg/webapp.yaml` or load your own config in `webapp/config.py`.
