InitBinder Bulk (Local Laptop Mode)
===================================

This repository includes the InitBinder web UI. For public GitHub usage, the recommended default is local laptop mode focused on:

- `http://127.0.0.1:8000/`
- If `8000` is already in use, `run_bulk_local.sh` automatically selects the next free local port
- Copy/paste command generation (no live cluster probing required)
- Optional advanced cluster execution handled by the user outside the app

Quickstart (macOS / Ubuntu)
---------------------------

Requires Python 3.10+.

1) Bootstrap a local environment:

```bash
./scripts/bootstrap_local.sh
```

This installs Python dependencies and Playwright Chromium used by vendor page parsing.

2) Optional preflight check:

```bash
./.venv/bin/python ./scripts/doctor_bulk_local.py
```

3) Run the server:

```bash
./scripts/run_bulk_local.sh
```

4) Open:

- `http://127.0.0.1:<selected_port>/` (defaults to `8000`, auto-shifts if occupied)

Catalog Input Preparation (Public Workflow)
-------------------------------------------

Recommended quick path (default for public users):

1) Use the prepared target table at:

- `targets_catalog/acrobio_plus_sino_biotin_merged.tsv`

2) Open `http://127.0.0.1:<selected_port>/` and paste TSV rows into **Input CSV / TSV**.

Optional full regeneration path (AcroBio + Sino + merge):

```bash
# Sino biotin catalog (raw + unique + manual)
python targets_catalog/webscraper/sino_biotin_pipeline.py

# AcroBio biotin catalog
python targets_catalog/webscraper/acrobio_biotin_pipeline.py \
  --url "https://www.acrobiosystems.com/search?keywords=biotinylated" \
  --mode playwright

# Merge two catalog TSVs into one bulk-ready file
python targets_catalog/webscraper/merge_target_catalogs.py \
  --input1 targets_catalog/webscraper/sino_biotinylated_unique.tsv \
  --input2 targets_catalog/webscraper/acrobio_biotinylated_unique.tsv \
  --output targets_catalog/acrobio_plus_sino_biotin_merged.tsv \
  --key uniprot
```

Notes on inputs/outputs:

- Inputs for generation are in `targets_catalog/webscraper/` (vendor scrape artifacts).
- Main merged output for bulk usage is `targets_catalog/acrobio_plus_sino_biotin_merged.tsv`.
- If you publish or redistribute scraped catalog data, verify vendor terms and applicable data-use policies first.

LLM API Requirements
--------------------

Some features require LLM APIs (for example, **Select epitopes (LLM)** in Bulk).

1) Open `Config` in the Bulk page.
2) Under `LLM`, set:
   - `OpenAI API key`
   - `OpenAI model` (optional, defaults to `gpt-4.1-mini`)
3) Save settings. These values persist to `cfg/webapp.local.yaml`.
4) Non-LLM actions (bulk preview, command generation, and most local/offline UI interactions) can still run without LLM keys.


Local-Safe Defaults
-------------------

The run script uses:

- `INITBINDER_UI_CONFIG=cfg/webapp.local.yaml` (if present) or `cfg/webapp.yaml`
- `INITBINDER_ALLOW_REMOTE=false`

This profile sets cluster mode to mock/offline defaults so command generation and local workflows do not depend on active SSH or remote root checks.


Remote Access Behavior
----------------------

By default, non-local clients are rejected (HTTP 403).  
To allow remote clients explicitly:

```bash
export INITBINDER_ALLOW_REMOTE=true
./scripts/run_bulk_local.sh
```

Use this only if you understand the network/security implications.


Command Generation (Bulk)
-------------------------

`Generate commands` and related command actions in Bulk are now config-driven and offline-friendly:

- Source = local config defaults + current UI overrides
- No cluster connection required to render commands
- Missing values use placeholders like `<cluster>` / `<remote_root>`


Configuration Precedence
------------------------

1. `cfg/webapp.yaml`
2. `cfg/webapp.local.yaml` (if present)
3. `INITBINDER_*` environment variables
4. If `INITBINDER_UI_CONFIG` is set, it overrides config path selection directly

For public/local mode, use `cfg/webapp.yaml` plus optional local overrides in `cfg/webapp.local.yaml`.


Notes
-----

- Bulk UI can run without PyMOL and without remote cluster access.
- LLM features require valid OpenAI credentials. For Bulk UI, set key/model in `Config` -> `LLM` (saved in `cfg/webapp.local.yaml`).
- Keep secrets and machine-specific overrides out of git (`cfg/webapp.local.yaml` is ignored).
- If any API key was committed or pasted into tracked files, rotate it immediately.
