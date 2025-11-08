# InitBinder Web UI

End-to-end FastAPI + Vanilla JS interface for orchestrating the InitBinder antibody design pipeline. The web layer wraps the existing CLI scripts (e.g. `manage_rfa.py`, `export_files.py`) so researchers can move from target setup to assessment review and PyMOL visualization without leaving the browser.

## Features

- **Target bootstrap** – queue `init-target`, `decide-scope`, and `prep-target` in a single click; live log streaming with job status.
- **Design orchestration** – push local inputs to the cluster, run `manage_rfa.py pipeline` remotely, submit RFdiffusion → ProteinMPNN → AlphaFold3 SLURM jobs, and auto-schedule `assess-rfa-all` once AF3 completes.
- **Results dashboard** – fetch AF3 rankings or BoltzGen metric CSVs, explore sortable tables, and inspect ipTM vs RMSD scatter plots with interactive highlighting.
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
  remote_root: /data/homezvol1/you/Projects/initbinder   # code/launcher location
  target_root: /pub/you/Projects/initbinder               # large design outputs
  pymol_path: pymol
  conda_activate: "source ~/.bashrc && conda activate takashi"
  assess_partition: standard             # optional CPU partition for post-processing
  assess_account: bio_lab                # optional SLURM account
  assess_time_minutes: 240               # walltime for assessment sbatch
  assess_mem_gb: 16
  assess_cpus: 4
  control_path: ~/.ssh/cm-initbinder     # optional explicit SSH control socket
  control_persist: 600                  # seconds to keep control master alive
  ensure_master: true
  mock: false
background_concurrency: 4
```

Environment overrides:

- `INITBINDER_PROJECT_ROOT` / `INITBINDER_ROOT`
- `INITBINDER_CLUSTER_ALIAS` / `INITBINDER_CLUSTER_HOST` / `INITBINDER_CLUSTER_USER`
- `INITBINDER_REMOTE_ROOT`
- `INITBINDER_TARGET_ROOT`
- `INITBINDER_CLUSTER_MOCK` (set to `true` for local dry-run)
- `INITBINDER_ASSESS_PARTITION`, `INITBINDER_ASSESS_ACCOUNT`, `INITBINDER_ASSESS_TIME_MINUTES`, `INITBINDER_ASSESS_MEM_GB`, `INITBINDER_ASSESS_CPUS`
- `INITBINDER_SSH_CONTROL_PATH`, `INITBINDER_SSH_CONTROL_PERSIST`, `INITBINDER_SSH_ENSURE_MASTER`
- `INITBINDER_CONDA_ACTIVATE`

### Cluster setup

- **Matching codebases** – the cluster path configured via `cluster.remote_root` must contain the *same git checkout* as the UI host (same branch/commit). The UI copies target inputs but assumes `manage_rfa.py` and generated tools already exist remotely.
- **SSH configuration** – create an entry in `~/.ssh/config` (e.g. `Host rfacluster`) that points to the SLURM login node. The web app shells out to `ssh`/`rsync`/`sbatch` using this alias. Example:

  ```sshconfig
  Host rfacluster
      HostName login.mycluster.edu
      User your_username
      Port 22
      ControlMaster auto
      ControlPath ~/.ssh/cm-initbinder-%r@%h:%p
      ControlPersist 10m
  ```
  Establish the control socket once per session (you will be prompted for your password or key the first time):

  ```bash
ssh hpc3.rcic.uci.edu -MNf
  ```

### How to completely remove SSH control master
If you need to reset the control master (e.g. after a network change), run:
```bash
cd ~/.ssh
rm cm-initbinder-*
```

  All subsequent `ssh`/`rsync` calls from the UI reuse this connection automatically.
- **Password-based logins** – if you cannot use keys, rely on the control master above so you only enter your password once per session. Tools like `sshpass` are discouraged; prefer SSH keys or control sockets.
- **Conda environment** – set `cluster.conda_activate` (e.g. `"source ~/.bashrc && conda activate takashi"`). Every remote Python command and assessment job prepends this snippet so the correct environment is active.
- **Environment modules** – ensure the remote environment can execute `python manage_rfa.py ...` and `sbatch`. The assessment helper creates jobs under `slurm_logs/` within the remote repository.
- **Automatic assessment** – after the RFdiffusion/MPNN/AF3 stage scripts are submitted, the UI schedules an additional `assess-rfa-all` sbatch job with dependencies on the final AF3 tasks. Configure the CPU partition/memory knobs above to match your cluster policies.

#### Assessment rescue behaviour

Occasionally a subset of stage-2 jobs (ProteinMPNN or AF3 batches) never exit
cleanly, leaving the dependent assessment job stuck in the `PD` state with
reasons such as `Dependency` or `DependencyNeverSatisfied`.  The design
workflow now launches a background *rescue monitor* whenever assessment was
requested:

- Every two minutes it snapshots SLURM status for the stage-2 job IDs.
- After three consecutive snapshots where **every** job remains pending solely
  for dependency-related reasons, it submits a fresh assessment job without any
  dependencies by calling `ClusterClient.submit_assessment(..., allow_empty_dependencies=True)`.
- All decisions (including the newly created sbatch ID) are logged to the job
  history so the UI timeline explains why assessment was resubmitted.

This keeps assessment results flowing even when some upstream stage-2 batches
fail or never launch, while preserving the normal dependency chain whenever
everything succeeds.

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
3. **Monitor runs** – Use SLURM tools externally or wait until results land. When ready, choose the result source (AF3 vs. BoltzGen), click the matching sync button (*Download assessments* or *Sync BoltzGen results*), then *Refresh* to pull the latest files.
4. **Analyze** – Sort/filter the table, click designs to highlight scatter points, and review metadata (epitope, PyMOL script path, CSV columns, etc.). When viewing BoltzGen metrics the scatter defaults to `design_to_target_iptm` vs. `filter_rmsd`.
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
| `GET` | `/api/targets/{pdb}/boltzgen/runs` | List local BoltzGen run/spec directories and CSV availability. |
| `POST` | `/api/targets/{pdb}/boltzgen/sync` | Rsync BoltzGen outputs from the cluster. |
| `GET` | `/api/targets/{pdb}/boltzgen/results` | Load BoltzGen `all_designs_metrics.csv` as a ranking payload. |
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
5. The UI schedules an `assess-rfa-all` job automatically (runs once the AF3 jobs succeed). After SLURM finishes, click *Sync from cluster* then *Refresh* to load the synthesized AF3 rankings.
6. Highlight top candidates via the scatter plot; inspect PyMOL script paths.
7. Set `Top N designs = 48`, codon host `yeast`, click *Generate exports* to produce ordering files.

---

For future enhancements and open issues, refer to `TODO.md`.
