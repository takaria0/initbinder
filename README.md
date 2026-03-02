InitBinder Bulk (Local Laptop Mode)
===================================

This repository includes the InitBinder web UI. For public GitHub usage, the recommended default is local laptop mode focused on:

- `http://127.0.0.1:8000/`
- If `8000` is already in use, `run_bulk_local.sh` automatically selects the next free local port
- LLM-first target discovery from a configured catalog
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

Bulk Workflow (LLM-First, Catalog-Backed)
-----------------------------------------

Recommended quick path:

1) In the Bulk page, open `Config`.
2) Set `Input TSV defaults` -> `Default input file path` to your catalog file.
3) Optionally enable `Auto-load default input`.
4) Use `LLM Target Discovery`:
   - choose/create a conversation,
   - describe the targets you want,
   - click `Suggest targets`.
5) Curate `Matched targets` with `Delete` / `Undelete`, and switch `View` (`LLM-picked` / `All targets`) plus `Catalog` filter (`Biotin only` / `All`) as needed.
6) Continue into `BoltzGen configs`, command generation, and binder review/export.

Notes:

- Default catalog used in this repo: `targets_catalog/acrobio_plus_sino_biotin_merged.tsv`.

After Selecting Targets
-----------------------

Practical next steps:

1) Run `Select epitopes (LLM)` and then `Rebuild configs`.
2) Validate epitopes using `Config` and `PyMOL`.
3) Open `Generate commands` (or per-row `Command`) to get `Cluster run steps`.
4) Execute generated commands manually in order, with local-vs-cluster boundaries:
   - `rsync`/`mkdir`/`ssh ...` monitor checks from local terminal
   - BoltzGen launch block in cluster shell
   - RFA launch block can run from local terminal as generated
5) Pull results, click `Refresh` in `Designed binders`, and review binder metrics (`ipTM`, `RMSD`, `Hotspot dist`, `ipSAE`).

For the full operator runbook, see `README_GUI.md`.

Catalog Refresh (Optional)
--------------------------

If you want to rebuild catalog files (AcroBio + Sino + merge):

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

LLM-backed features include both target discovery and epitope actions, including:

- `LLM Target Discovery` (`Suggest targets`)
- `Select epitopes (LLM)` actions in Bulk tables/modals

Setup:

1) Open `Config` in the Bulk page.
2) Under `LLM`, set:
   - `OpenAI API key`
   - `OpenAI model` (optional, defaults to `gpt-4.1-mini`)
3) Save settings. These values persist to `cfg/webapp.local.yaml`.
4) Non-LLM actions (preview, config generation, command generation, and local table review) can still run without LLM keys.

Export Selected Binders (Current Behavior)
------------------------------------------

`Export selected binders` now produces more than a plain binder CSV:

- Selects top binders per `PDB:epitope` using combined `ipTM + RMSD` ranking.
- Adds yeast-codon DNA plus deterministic constrained adapters.
  - Binders sharing the same `PDB + hotspot` use the same adapter seed/barcode.
- Runs BsaI checks on `adapter + binder` combined sequence.
  - Automatically retries codon optimization to reduce extra BsaI sites.
  - Unresolved rows are still exported and flagged with `bsai_site_check_ok=False`.
- Exports scatter artifacts per engine (`ipTM` vs `RMSD`) with antigen:epitope color grouping.

Prompt Examples
---------------

For ready-to-use LLM discovery prompt ideas, see:

- `example_prompts.md`

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

`Generate commands` and related command actions in Bulk are config-driven and offline-friendly:

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
- If any API key was committed or entered into tracked files, rotate it immediately.
