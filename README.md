InitBinder Bulk (Local Laptop Mode)
===================================

This repository includes the InitBinder web UI. For public GitHub usage, the recommended default is local laptop mode focused on:

- `http://127.0.0.1:8000/bulk`
- Copy/paste command generation (no live cluster probing required)
- Optional advanced cluster execution handled by the user outside the app


Quickstart (macOS / Ubuntu)
---------------------------

1) Bootstrap a local environment:

```bash
./scripts/bootstrap_local.sh
```

2) Optional preflight check:

```bash
./scripts/doctor_bulk_local.py
```

3) Run the server:

```bash
./scripts/run_bulk_local.sh
```

4) Open:

- `http://127.0.0.1:8000/bulk`


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
- LLM features require valid API keys. Copy `cfg/env.sample.py` to `cfg/env.py` and set keys.
- Keep secrets and machine-specific overrides out of git (`cfg/env.py`, `cfg/webapp.local.yaml` are already ignored).
