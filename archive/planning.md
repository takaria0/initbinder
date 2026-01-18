# UI Revamp Planning

## Requested UI Improvements (Apr 2025)
- Introduce a debug/advanced toggle that reveals `Force Refresh`, `Run Decide-Scope`, and `Run Prep-Target`; hide these controls during normal use.
- Extend target setup form with an explicit `Number of Epitopes` input tied to downstream pipeline configuration.
- In Sequence Alignment view, display both termini of unmatched regions alongside the mismatched length for clearer context.
- Within Design Run Configuration, remove the `Binder Chain ID` input from the primary UI (move to debug/advanced if still needed).
- In Assessment Results:
  - Hide the run-label input and rename `Limit` to `Row Limit`.
  - Remove visible `Reload Local Data` and `Download from Cluster` buttons; instead poll automatically every few seconds (these controls can live under the debug toggle if required).
  - Move the `PyMOL Top 96` button to the top control row of the TSV viewer (before table header) and add a sibling `Export` button that opens a modal for configuring exports.
- Add persistent tracking of PDB ID + Antigen URL pairs (with user-defined names) so users can quickly re-select past combinations; selecting a saved pair should auto-populate the relevant fields throughout the app.
- Adjust gallery handling so `assess_all` (or equivalent) writes gallery outputs into the same directory that hosts `af3_rankings.tsv`; update the UI so the PyMOL gallery launcher targets the colocated folder.

### Design Notes
- Introduce a `debugMode` flag in frontend state driven by a toggle; hide advanced controls with a shared CSS class when disabled.
- Add a new backend store (`TargetPreferences`) that keeps named PDB/antigen pairs on disk (`cache/ui_state/targets.json`) and expose `/api/targets/prefs` endpoints for list/create/update/delete.
- Extend `TargetInitRequest` schema with `num_epitopes` and persist to `target.yaml` through existing workflow helpers.
- Sequence alignment API should return mismatched flank residues (`left_context`/`right_context`) in addition to length; update renderer to surface them.
- Rework rankings polling so the UI refreshes via interval when the table is visible; manual reload/sync buttons move into debug panel.
- `designs.run_design` should treat binder chain as optional; only use value when provided.
- Create modal infrastructure (simple hidden div) to support TSV export configuration triggered by the new top-level Export button.
- Attach PDB/antigen selections to design + assessment requests so dependent forms auto-populate when a saved set is chosen.
- Update assessment workflow so gallery artifacts are copied into the same leaf directory as `af3_rankings.tsv`; adjust PyMOL launcher to reference that path.


## Goals
- Provide a single-page web UI where a scientist can feed a PDB ID and vendor antigen URL, then drive the entire InitBinder pipeline without touching the CLI.
- Expose design artifacts visually (sequence alignment, plots, tables) and operationally (cluster batch submissions, PyMOL launches, export bundles).
- Keep the implementation modular so we can reuse existing Python pipeline scripts while presenting a modern HTML/CSS/JS experience.

## High-Level Workflow
1. **Target Intake**
   - User supplies PDB ID + antigen URL.
   - Backend invokes `manage_rfa.py init-target` → writes/updates `targets/<PDB>/target.yaml`.
   - Run `manage_rfa.py decide-scope` and `manage_rfa.py prep-target` automatically (with status streaming).
   - Cache the resulting `target.yaml`, epitope masks, hotspot variants, etc. for downstream steps.

2. **Visualization & QA**
   - Align antigen sequence (from vendor URL or `target.yaml` accession) vs. PDB chains using `sequence_alignment.py`.
   - Render alignment and epitope annotations inside the browser (selectable mismatch view).
   - Offer “Launch PyMOL” button to open hotspot bundle locally via `scripts/pymol_utils.py`.

3. **Design Execution (Cluster)**
   - User chooses number of designs per engine (RFdiffusion, ProteinMPNN, AF3) and run label.
   - Backend packages target folder, rsyncs/ssh’s to HPC, generates SLURM scripts (`make_rfa_*` helpers), optionally submits jobs, tracks job IDs + log paths.
   - Poll job status (via `squeue`) and show progress.

4. **Result Aggregation**
   - Once designs finish, pull back new `af3_rankings.tsv`, PyMOL bundles, and logs.
   - Run `plot_rankings.py` to get scatter plot JSON/PNG; serve data for front-end interactive scatter (IPTM vs RMSD, clickable points revealing binder metadata, direct PyMOL launch for the specific binder).
   - Render `af3_rankings.tsv` in sortable/filterable table. Provide `order_by` toggles.

5. **Exports**
   - Trigger `export_files.py` from UI with user-selected parameters (top N, codon options). Offer downloads (FASTA, CSV, Excel).
   - For binder selection, allow multi-select + “Launch PyMOL (Top N)” button to stream PML to local PyMOL.

## Proposed Architecture

### Backend (Python + FastAPI)
- New package under `webapp/` (since historical files are gone) with structure:
  ```
  webapp/
    __init__.py
    config.py              # Cluster + path settings (read from env/local YAML)
    main.py                # FastAPI app entrypoint
    hpc.py                 # SSH/rsync helpers, SLURM job orchestration
    pipeline.py            # Wrappers around manage_rfa / make_rfa_* scripts
    alignment.py           # Sequence alignment utilities → JSON for UI
    result_collectors.py   # Ranking table, plot generation, export invocations
    pymol.py               # Launch + bundle generation helpers
    models.py              # Pydantic schemas for requests/responses
    job_store.py           # In-memory + on-disk job tracking (status, logs)
    workflows.py           # Higher-level orchestrations (init + prep pipeline, design run, assess, export)
    static/                # Bundled JS/CSS assets served by FastAPI
    templates/             # `index.html`
  ```
- Use `uvicorn webapp.main:app --reload` during development; production can reverse-proxy.
- Long-running tasks handled via `BackgroundTasks` or `async` workers; statuses kept in `job_store` for polling.
- HPC credentials taken from `~/.ssh/config` alias + environment variables; support overrides via UI form or config file.

### Frontend (Vanilla JS + HTMX-style components)
- Single-page app delivered from `templates/index.html`, minimal build tooling.
- Modules for:
  - Target form & status feed (progress log stream via EventSource/`/stream/{job_id}` endpoint).
  - Sequence alignment canvas (custom JS/Canvas or simple HTML table with color-coded spans).
  - PyMOL controls (`Launch Hotspots`, `Launch Top N Binders`).
  - Design configuration panel (sliders/inputs for counts; checkboxes to include AF3/MPNN/RFdiff).
  - Results dashboard: sortable table (use `Tabulator` or lightweight custom sorter), scatter plot (use Plotly.js for interactive points), downloads section.
- Keep CSS modular (`static/css/app.css`); adopt CSS grid/flex for layout.

### HPC Integration
- Parse SLURM submission info returned by `make_rfa_*` helpers to capture job scripts + names.
- Provide config for remote root (e.g., `/pub/inagakit/Projects/initbinder`). Tools:
  - `rsync` for directory sync (target push + result pull).
  - `ssh` to run remote commands: module loads, `sbatch`, `squeue`, `sacct`.
- Maintain mapping `job_id -> {stage: str, state: pending/running/done, logs: list[str], remote_paths: ...}`.
- Offer UI button to “Download Now” once `assess` stage completes (fires rsync).

## Feature Checklist & Implementation Notes

- [ ] Target initialization form → triggers `POST /targets/{pdb_id}` (background job) performing init → decide-scope → prep pipeline.
- [ ] Sequence alignment endpoint → returns JSON with aligned rows, mismatch annotations, epitope overlays.
- [ ] PyMOL hotspot launch → local exec of existing bundle helper, optionally direct `pymol` invocation if installed.
- [ ] Design run form → `POST /designs` with payload (target, arms, counts, binder chain, run label, cluster opts). Backend:
      1. Ensure local prep complete.
      2. Sync to remote HPC path.
      3. Invoke `make_rfa_rfdiffusion`, `make_rfa_proteinmpnn`, `make_rfa_af3`, `assess_rfa_design` or `manage_rfa pipeline` as appropriate (generate scripts & optionally submit).
      4. Track SLURM job IDs returned by `sbatch` (parse output).
- [ ] Status polling route `/jobs/{job_id}` returning human-friendly updates (progress text, HPC job states, next actions).
- [ ] Result ingestion `GET /targets/{pdb_id}/rankings` reading latest TSV (first local, else remote fetch) and delivering JSON.
- [ ] Scatter plot data transformation (call `plot_rankings.py` with `--json_out` patch if necessary, otherwise parse TSV and compute on backend, deliver to UI; still run script to save static plot for download).
- [ ] Data table sorting/filtering in JS (bring your own logic; `Intl.Collator` for typed comparisons; maintain order_by state).
- [ ] Bindings to run `export_files.py` with UI-specified parameters; store outputs under `cache/ui_exports/<timestamp>/` and return file download links.
- [ ] Optional: integrate websockets/EventSource for log streaming from backend to UI for long operations.

## Dependencies & Tooling
- Python: FastAPI, Uvicorn, `paramiko` (or `asyncssh`) for cluster SSH, `aiofiles` for streaming logs, `pydantic` already included with FastAPI.
- Frontend: Plotly.js (self-hosted via `/static/js/plotly.min.js`), maybe `Sortable` for table columns; keep dependencies vendored (no CDN due to offline clusters).
- For background execution, rely on `asyncio.create_task` + `job_store` (thread-safe) or use `anyio` tasks built into FastAPI.
- Ensure `requirements.txt` captures new packages; add `webapp/__init__.py` to mark package.

## Risks & Mitigations
- **ROOT path mismatch** (`utils.ROOT` hardcoded to `/pub/...`): wrap helper to temporarily set env var or patch utils to detect local vs HPC root per config.
- **Long-running tasks**: guard with timeouts, send periodic heartbeats to UI.
- **PyMOL invocation**: run on local workstation; ensure backend checks for binary presence and provides clear error messages.
- **Cluster credentials**: avoid storing passwords; rely on SSH keys configured in `~/.ssh/config`. Provide UI to choose host alias (default `rfacluster`).
- **Plotting**: `plot_rankings.py` currently writes files; extend to emit JSON without heavy rewrites (maybe wrap in backend parser).

## Implementation Milestones
1. **Backend skeleton**: FastAPI app, job store, config loader, CLI wrappers for init/decide/prep.
2. **Frontend shell**: HTML layout, fetch hooks, progress UI, alignment widget.
3. **Cluster module**: SSH/rsync wrappers, job submission/poll API.
4. **Design orchestration**: end-to-end run (RFdiff → MPNN → AF3 → Assess) with status updates.
5. **Results dashboard**: load rankings, render table/chart, provide downloads.
6. **Polish**: PyMOL integration, export workflows, error handling, documentation.

## Testing Strategy
- Local integration test endpoints using sample target (e.g., 6M17) with mocked HPC interactions (flag in config to simulate cluster via local subprocess).
- Unit-test alignment JSON + ranking parser with fixtures from `targets/6M17`.
- Manual UI smoke test via `uvicorn` and mock cluster to validate flows.
- Provide `make webapp-dev` convenience script (optional) for running server with auto reload.

## Open Questions
- Should design runs use `manage_rfa.py pipeline` (single command) or the per-stage scripts? (Plan: start with `pipeline` for coarse control, add advanced toggles later.)
- For scatter plot interactivity, is Plotly acceptable or should we craft custom D3? (Default: Plotly for speed.)
- Export settings default values (GC target, codon host) – expose minimal necessary knobs initially.
