# InitBinder Web UI

End-to-end FastAPI + Vanilla JS interface for orchestrating the InitBinder antibody design pipeline. The web layer wraps the existing CLI scripts (e.g. `manage_rfa.py`, `export_files.py`) so researchers can move from target setup to assessment review and PyMOL visualization without leaving the browser.

## Features

- **Target bootstrap** – queue `init-target`, `decide-scope`, and `prep-target` in a single click; live log streaming with job status.
- **Design orchestration** – submit multi-arm RFdiffusion → ProteinMPNN → AlphaFold3 batches to the cluster using the same heuristics as the CLI launcher.
- **Results dashboard** – fetch latest AF3 rankings, explore sortable tables, and inspect ipTM vs RMSD scatter plots with interactive highlighting.
- **PyMOL automation** – generate hotspot bundles and top-binder galleries (local or remote `pymol-remote` modes supported).
- **Export utilities** – drive `export_files.py` directly from the UI to create FASTA/CSV/Excel order packets.
- **Cluster sync** – rsync assessment outputs back from the HPC and rehydrate the UI for post-run analysis.

## Architecture Overview

```
Browser (HTML + JS)
  │
  ├─ REST API (FastAPI, `webapp/main.py`)
  │    ├─ Job store (`job_store.py`) for status & log aggregation
  │    ├─ Workflows (`workflows.py`) dispatch background threads
  │    ├─ Pipeline wrappers (`pipeline.py`, `designs.py`, `exporter.py`)
  │    ├─ Data loaders (`alignment.py`, `result_collectors.py`)
  │    └─ PyMOL helpers (`pymol.py`, `scripts/pymol_utils.py`)
  │
  └─ HPC client (`hpc.py`) handles rsync/ssh with mock or real cluster
```

### Backend packages

| Module | Responsibility |
| --- | --- |
| `config.py` | Load config from YAML/env (paths, cluster settings, concurrency). |
| `job_store.py` | Thread-safe in-memory log/status persistence with optional JSON snapshots. |
| `pipeline.py` | Thin wrapper around `manage_rfa.py` subcommands. |
| `designs.py` | Writes SLURM launchers, syncs project data, and submits to cluster via SSH. |
| `hpc.py` | rsync/ssh utilities (mock mode stores data under `.cluster-mock/`). |
| `result_collectors.py` | Parse AF3 ranking TSVs for UI consumption. |
| `pymol.py` | Generate and optionally launch PyMOL bundles (hotspots, top-N binders). |
| `exporter.py` | Stream `export_files.py` output into the job store. |
| `static/` | Pure JS/CSS UI (no build tooling). |

### Frontend structure

- `templates/index.html` serves a single SPA-style page
- `static/js/app.js` handles fetch requests, state management, job polling, scatter plot rendering, and table interactivity
- `static/css/app.css` provides a lightweight design system (no external frameworks)

## Installation

1. **Install Python deps (UI server + collectors):**
   ```bash
   pip install -r requirements-webapp.txt
   ```
   The UI relies on packages already used by the CLI (`pandas`, `numpy`, `dnachisel`, etc.)—ensure the main InitBinder environment exists.

2. **Optional PyMOL integration:**
   - Local mode: install PyMOL on the workstation and ensure `pymol` is on `PATH`.
   - Remote mode: `pip install pymol-remote` on both the cluster and the local machine, then set `RFA_PYMOL_MODE=remote`.

3. **Cluster access:** configure passwordless SSH to your SLURM login node. The backend uses `ssh`, `rsync`, and `sbatch` via the alias described in `webapp/config.py` / `cfg/webapp.yaml`.

## Configuration

Create `cfg/webapp.yaml` (or set environment variables) to match your environment. Example:

```yaml
paths:
  project_root: /Users/you/Projects/initbinder
cluster:
  ssh_config_alias: rfacluster            # matches ~/.ssh/config host entry
  remote_root: /pub/you/Projects/initbinder
  pymol_path: pymol
  mock: false
background_concurrency: 4
```

Environment overrides:

- `INITBINDER_PROJECT_ROOT` / `INITBINDER_ROOT`
- `INITBINDER_CLUSTER_ALIAS` / `INITBINDER_CLUSTER_HOST` / `INITBINDER_CLUSTER_USER`
- `INITBINDER_REMOTE_ROOT`
- `INITBINDER_CLUSTER_MOCK` (set to `true` for local dry-run)

## Running the Server

```bash
# From repository root
uvicorn webapp.main:app --reload --port 8000
```

Visit <http://localhost:8000/> to load the UI. The landing page serves static assets directly from `templates/` and `static/`.

### Development helpers

```bash
# Lint the backend (optional)
flake8 webapp

# Regenerate type hints (optional)
pyright webapp
```

## Typical Workflow

1. **Target setup** – Enter PDB ID + antigen URL, click *Queue Pipeline*. Watch logs for `init-target`, `decide-scope`, and `prep-target` completion. Once done, PyMOL hotspots become available.
2. **Design submission** – Paste epitope arms (`<Epitope>@<Variant>`), adjust counts, set run label, and submit. The backend writes/patches the launcher, syncs to the cluster, and executes via SSH.
3. **Monitor runs** – Use SLURM tools externally or wait until results land. When ready, *Sync from cluster* and *Refresh* to pull the latest AF3 ranking TSV.
4. **Analyze** – Sort/filter the rankings table, click designs to highlight scatter points, and review metadata (epitope, PyMOL script path, etc.).
5. **PyMOL review** – Launch hotspot bundles or top binders (96 by default, adjustable via Export Top N).
6. **Export** – Set FASTA/CSV options (codon host, GC target, prefix/suffix) and generate full ordering packets.

## API Endpoints (brief)

| Method | Path | Description |
| --- | --- | --- |
| `POST` | `/api/targets/init` | Run `init-target` → `decide-scope` → `prep-target`. |
| `POST` | `/api/designs/run` | Generate scripts, sync to cluster, run SLURM pipeline. |
| `POST` | `/api/exports` | Execute `export_files.py` with provided options. |
| `GET` | `/api/jobs/{id}` | Retrieve latest job logs/status. |
| `GET` | `/api/targets/{pdb}/alignment` | Sequence alignment summary from `target.yaml`. |
| `GET` | `/api/targets/{pdb}/rankings` | Load parsed AF3 rankings TSV. |
| `POST` | `/api/targets/{pdb}/pymol/hotspots` | Build (and optionally launch) hotspot bundle. |
| `POST` | `/api/targets/{pdb}/pymol/top-binders` | Build Multi-binder PyMOL session from rankings. |
| `POST` | `/api/targets/{pdb}/sync` | Rsync assessments from cluster to local workspace. |

All long operations respond immediately with a job ID; the UI polls `/api/jobs/{id}` to surface progress.

## Mock Mode

Set `INITBINDER_CLUSTER_MOCK=true` (or `mock: true` in the config) to test without a cluster. The HPC client mirrors rsync targets under `<project>/.cluster-mock`, while SSH calls become no-ops.

## Troubleshooting

- Ensure `INITBINDER_ROOT` or `paths.project_root` match the repository checkout; otherwise the backend may not find `manage_rfa.py`.
- If PyMOL buttons fail, confirm `scripts/pymol_utils.py` imports successfully and PyMOL is on `PATH`.
- For missing rankings, verify `targets/<PDB>/designs/_assessments/<run>/af3_rankings.tsv` exists locally or run the sync step.
- For cluster submission errors, inspect job logs for `[sbatch]` outputs (the log stream includes SLURM job IDs).

## Example Session

1. `uvicorn webapp.main:app --reload`
2. Navigate to the UI → `PDB ID = 6M17`, enter antigen URL, queue pipeline.
3. After prep completes, click *Launch hotspots in PyMOL* to view hotspots locally.
4. Paste arm list (e.g. `RBM Core@A`, `RBM Core@B`), set `Run label = demo`, submit design pipeline.
5. Once SLURM jobs finish, click *Sync from cluster* then *Refresh* to load AF3 rankings.
6. Highlight top candidates via the scatter plot; inspect PyMOL script paths.
7. Set `Top N designs = 48`, codon host `yeast`, click *Generate exports* to produce ordering files.

---

For future enhancements and open issues, refer to `TODO.md`.

