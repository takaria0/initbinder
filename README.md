InitBinder
==========

[Demo Video](https://www.youtube.com/watch?v=ILdtOMNKADo)

InitBinder is a local web application for catalog-guided antigen selection, epitope curation, cluster command generation, binder review, and synthesis ready binder export in antibody design workflows. The public entrypoint in this repository is served at `http://127.0.0.1:8000/`.

Users can discover purchasable targets from a configured catalog, select epitopes with LLM-assisted or manual workflows, generate BoltzGen configs, generate copy/paste cluster commands, and export ranked binders with DNA/adaptor annotations. For the detailed operator runbook, see `README_GUI.md`.

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

This creates `.venv`, installs the runtime dependencies from `requirements-webapp.txt`, and installs Playwright Chromium used by the local tooling (Note: This does not have to be venv, if you know what you're doing, go ahead and use your preferred virtual environment setup).

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
3. In `Cluster`, set `Cluster user` and `Cluster host` if you want to use the `Generate commands` feature to prepare remote execution commands. This does not require a live SSH connection just to render the commands, but you will need cluster access to execute them later. You also need to install BoltzGen and this application on the cluster if you want to run the generated commands remotely.
4. In `LLM`, add an `OpenAI API key` if you want to use `Suggest targets` or `Select epitopes (LLM)`.
5. In `LLM Target Discovery`, create or select a conversation, describe the targets you want, and click `Suggest targets`.
6. Curate `Matched targets` with `Delete` / `Undelete`, then confirm the retained rows in `Selected targets`.
7. In `BoltzGen configs`, set `Epitope design num`, `Max runtime (hours)`, and `Crop radius`.
8. Run `Select epitopes (LLM)` when needed, or add manual epitopes from the per-target `Select epitopes` modal.
9. Wait for the boltzgen config to be generated, or click `Rebuild configs` to get boltzgen configs.
10. Use `Generate commands` or per-row `Command` to open `Cluster run steps`.
11. Execute the generated command blocks manually, sync results back, click `Refresh` in `Designed binders`, and export the selected binders when ready.


LLM and cluster requirements
----------------------------

LLM-backed actions include:

- `LLM Target Discovery` -> `Suggest targets`
- `Select epitopes (LLM)` in the default workflow

These features require a valid OpenAI API key configured in `Config` -> `LLM`. Non-LLM operations, including config regeneration, command rendering, local review, and most export flows, can still run without LLM credentials.

Cluster access is optional for command generation but required for actual remote design execution. The UI generates commands from the saved config and current UI state; it does not need a live SSH connection just to render the command blocks.


Limitations
-----------

- The public workflow is local-first and assumes a trusted single-user workstation.
- Command execution is manual; the application prepares commands but does not fully orchestrate remote jobs end-to-end.
- Actual design execution depends on external cluster access, installed design software, and user-managed credentials.
- LLM suggestions are constrained by the configured catalog and by the quality of the prompt.

