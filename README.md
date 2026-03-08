InitBinder
==========

InitBinder is a local-first web application for catalog-guided antigen selection, epitope curation, cluster command generation, and binder review/export in antibody design workflows. The public entrypoint in this repository is served at `http://127.0.0.1:8000/`, with local-safe defaults that avoid live cluster probing and reject remote clients unless explicitly enabled.

The current implementation emphasizes interactive target triage: users can discover purchasable targets from a configured catalog, select epitopes with LLM-assisted or manual workflows, regenerate BoltzGen inputs, generate copy/paste cluster commands, and export ranked binders with DNA/adaptor annotations. For the detailed operator runbook, see `README_GUI.md`.

Availability and implementation
-------------------------------

- Language: Python 3.10+
- Application shape: FastAPI backend with a browser-based interface
- License: MIT (`LICENSE`)
- Default runtime mode: localhost-only, single-user, local laptop workflow
- Default catalog in this repository: `targets_catalog/acrobio_plus_sino_biotin_merged.tsv`

Installation
------------

Bootstrap the runtime environment:

```bash
./scripts/bootstrap_local.sh
```

This creates `.venv`, installs the runtime dependencies from `requirements-webapp.txt`, and installs Playwright Chromium used by the local tooling.

Optional checks and add-on installs:

```bash
# Preflight check for the runtime
./.venv/bin/python ./scripts/doctor_bulk_local.py

# Test dependencies
./.venv/bin/python -m pip install -r requirements-dev.txt

# Optional utility dependencies for catalog refresh / diversity extras
./.venv/bin/python -m pip install -r requirements-optional.txt
```

Start the server:

```bash
./scripts/run_bulk_local.sh
```

Then open `http://127.0.0.1:<selected_port>/`. The launcher defaults to port `8000` and automatically shifts to the next free local port if needed.

Quickstart
----------

1. Open `Config`.
2. In `Input TSV defaults`, set `Default input file path` to your catalog TSV/CSV and optionally enable `Auto-load default input`. To follow the paper's default reproduction path in this repository, use `targets_catalog/acrobio_plus_sino_biotin_merged.tsv`; it already includes the required entries, including the Sino Biological human TIM-3 antigen (`10390-H08H-B`).
3. In `LLM`, add an `OpenAI API key` if you want to use `Suggest targets` or `Select epitopes (LLM)`.
4. In `LLM Target Discovery`, create or select a conversation, describe the targets you want, and click `Suggest targets`.
5. Curate `Matched targets` with `Delete` / `Undelete`, then confirm the retained rows in `Selected targets`.
6. In `BoltzGen configs`, set `Epitope design num`, `Max runtime (hours)`, and `Crop radius`.
7. Run `Select epitopes (LLM)` when needed, or add manual epitopes from the per-target `Select epitopes` modal.
8. Click `Rebuild configs`.
9. Use `Generate commands` or per-row `Command` to open `Cluster run steps`.
10. Execute the generated command blocks manually, sync results back, click `Refresh` in `Designed binders`, and export the selected binders when ready.

Inputs and outputs
------------------

Typical inputs:

- Catalog TSV/CSV containing `chosen_pdb` or `pdb_id`
- Optional vendor metadata such as `antigen_url`, `vendor_accession`, and `vendor_range`
- Optional OpenAI credentials for LLM-assisted actions
- Optional cluster settings for actual remote execution

Primary outputs:

- Regenerated BoltzGen configuration files for active targets/epitopes
- Copy/paste cluster execution commands
- Epitope and diversity reports (`HTML` / `CSV`)
- Binder review tables with ranking metrics such as `ipTM`, `RMSD`, `Hotspot dist`, and `ipSAE`
- `Export selected binders` artifacts, including ranked binder CSV rows, codon-optimized DNA, deterministic adapter annotations, and per-engine scatter plots

LLM and cluster requirements
----------------------------

LLM-backed actions include:

- `LLM Target Discovery` -> `Suggest targets`
- `Select epitopes (LLM)` in the default workflow

These features require a valid OpenAI API key configured in `Config` -> `LLM`. Non-LLM operations, including config regeneration, command rendering, local review, and most export flows, can still run without LLM credentials.

Cluster access is optional for command generation but required for actual remote design execution. The UI generates commands from the saved config and current UI state; it does not need a live SSH connection just to render the command blocks.

Reproducibility and testing
---------------------------

Runtime installation is defined by `requirements-webapp.txt`. Test execution is defined by `requirements-dev.txt`. Optional catalog-refresh and plotting extras are separated into `requirements-optional.txt` so the default install remains lightweight.

Recommended local verification commands:

```bash
# Runtime import / syntax sanity
python3 -m py_compile webapp/*.py lib/scripts/*.py *.py targets_catalog/webscraper/*.py

# Test suite
./.venv/bin/python -m pytest -q
```

The runtime launcher uses:

- `INITBINDER_UI_CONFIG=cfg/webapp.local.yaml` if present, otherwise `cfg/webapp.yaml`
- `INITBINDER_ALLOW_REMOTE=false` by default

Remote access behavior
----------------------

By default, non-local clients receive HTTP `403`. To allow remote clients explicitly:

```bash
export INITBINDER_ALLOW_REMOTE=true
./scripts/run_bulk_local.sh
```

Use this only if you are prepared to provide your own network isolation, authentication, and TLS termination.

Limitations
-----------

- The public workflow is local-first and assumes a trusted single-user workstation.
- Command execution is manual; the application prepares commands but does not fully orchestrate remote jobs end-to-end.
- Actual design execution depends on external cluster access, installed design software, and user-managed credentials.
- LLM suggestions are constrained by the configured catalog and by the quality of the prompt.
- If you rebuild or redistribute scraped catalog content, verify vendor terms and local data-use requirements first.

Security and data use
---------------------

- Keep secrets and machine-specific overrides out of git. `cfg/webapp.local.yaml` is ignored for that reason.
- If any API key was previously committed or stored in tracked files, rotate it immediately.
- See `SECURITY.md` for the local-mode security posture and the expectations if remote access is enabled.
