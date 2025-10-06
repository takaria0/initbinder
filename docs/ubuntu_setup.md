# InitBinder Ubuntu Setup Guide

This document captures the full end-to-end steps to bootstrap the InitBinder
antibody design environment on a fresh Ubuntu workstation so the CLI tools and
web UI both function and can drive the remote HPC workflow.

## 1. System prerequisites

1. Update the package index and install core build tooling, Git, and Python 3.10
   with its virtual environment helpers:
   ```bash
   sudo apt-get update
   sudo apt-get install -y git build-essential python3.10 python3.10-venv python3-pip rsync openssh-client
   ```
2. (Optional) Install PyMOL locally if you plan to launch hotspot or top-binder
   bundles directly from the UI. On Ubuntu this is typically:
   ```bash
   sudo apt-get install -y pymol
   ```
   You can skip this on headless servers and rely on exported bundles instead.

## 2. Clone the repository

```bash
git clone git@github.com:allen-cell-immune/initbinder.git
cd initbinder
```

If you prefer HTTPS:
```bash
git clone https://github.com/allen-cell-immune/initbinder.git
```

> The CLI scripts reference files relative to the repository root, so keep the
> checkout path stable or adjust `ROOT` inside `utils.py` accordingly.

## 3. Python virtual environment

1. Create and activate a dedicated virtual environment (Python ≥3.10):
   ```bash
   python3.10 -m venv .venv
   source .venv/bin/activate
   python -m pip install --upgrade pip
   ```
2. Record the activation command for later (e.g. add `source /path/to/.venv/bin/activate`
   to your shell profile).

## 4. Install Python dependencies

### 4.1 CLI & discovery pipeline

Install the libraries needed for discovery, structure preparation, and the
assessment stages. These are enumerated in the root `README.md`.

```bash
pip install google-generativeai playwright bs4 requests pyyaml tenacity httpx fake_useragent \
            biopython freesasa jsonschema numpy pandas dnachisel
python -m playwright install chromium
```

- Ensure `env.py` exports `GOOGLE_API_KEY` for LLM-powered discovery (`target_generation.py`).
- `freesasa` installs compiled extensions; if the wheel build fails, install the
  system headers via `sudo apt-get install -y libfreesasa-dev` and retry.

### 4.2 Web UI backend

Install the FastAPI server requirements that wrap the CLI:
```bash
pip install -r requirements-webapp.txt
```

Optional extras:
- `pip install pymol-remote` on both local and remote machines when using
  `RFA_PYMOL_MODE=remote`.
- `pip install uvloop` to accelerate the ASGI server.

## 5. Environment configuration

### 5.1 Repository settings

1. Copy and edit `cfg/webapp.yaml` to match your workstation and cluster paths.
   Critical fields (see `webapp/README.md`) include:
   ```yaml
   paths:
     project_root: /home/you/Projects/initbinder
   cluster:
     ssh_config_alias: rfacluster
     remote_root: /data/homezvol1/you/Projects/initbinder
     target_root: /pub/you/Projects/initbinder
     conda_activate: "source ~/.bashrc && conda activate takashi"
     assess_partition: standard
     assess_account: bio_lab
     assess_time_minutes: 240
     assess_mem_gb: 16
     assess_cpus: 4
     control_path: ~/.ssh/cm-initbinder
     control_persist: 600
     ensure_master: true
     mock: false
   background_concurrency: 4
   ```
2. If your local checkout lives elsewhere, edit `ROOT` and related constants in
   `utils.py` so generated scripts and rsync paths resolve correctly.

### 5.2 Environment variables

Export these from your shell profile or `.env` file:
```bash
export GOOGLE_API_KEY="<your_key>"
export INITBINDER_PROJECT_ROOT="/home/you/Projects/initbinder"
export INITBINDER_CLUSTER_ALIAS="rfacluster"
export INITBINDER_REMOTE_ROOT="/data/homezvol1/you/Projects/initbinder"
export INITBINDER_TARGET_ROOT="/pub/you/Projects/initbinder"
export INITBINDER_CONDA_ACTIVATE="source ~/.bashrc && conda activate takashi"
```
Set `INITBINDER_CLUSTER_MOCK=true` to dry-run without a cluster.

## 6. SSH and cluster preparation

1. **Generate a key pair (if you do not already have one):**
   ```bash
   ssh-keygen -t ed25519 -C "you@initbinder"
   eval "$(ssh-agent -s)"
   ssh-add ~/.ssh/id_ed25519
   ```
2. **Authorize the key on the cluster login node:**
   ```bash
   ssh-copy-id your_username@login.mycluster.edu
   ```
3. **Create a host alias in `~/.ssh/config` for the web app to reference:**
   ```sshconfig
   Host rfacluster
       HostName login.mycluster.edu
       User your_username
       Port 22
       IdentityFile ~/.ssh/id_ed25519
       ControlMaster auto
       ControlPath ~/.ssh/cm-initbinder-%r@%h:%p
       ControlPersist 10m
   ```
4. **Prime the control master socket before launching the UI:**
   ```bash
   ssh -MNf rfacluster
   ```
   You should only be prompted for your password or key passphrase once per
   session.
5. Ensure `rsync`, `ssh`, and `sbatch` are available on the cluster, and that the
   remote path configured in `cfg/webapp.yaml` contains a matching Git checkout
   (`git clone` + `git pull` on the cluster as needed).
6. Confirm the remote Conda environment (specified by `cluster.conda_activate`)
   can run `python manage_rfa.py ...` and that SLURM partitions referenced in
   the config exist.

## 7. Typical CLI workflow

Run the orchestrator commands from the repository root after activating your
virtual environment:
```bash
# Discovery
python manage_rfa.py target-generation --instruction "top proteins for antibody therapeutics" \
  --max_targets 100 --species human --prefer_tags biotin --no_browser_popup

# Target setup
python manage_rfa.py init-target 6M17 --antigen_url "https://www.sinobiological.com/..."
python manage_rfa.py decide-scope 6M17
python manage_rfa.py prep-target 6M17 --sasa_cutoff 10.0

# Design pipeline (edit arms to match prep output)
python manage_rfa.py pipeline 6M17 \
  --arm "RBM Core@A" --arm "RBM Core@B" --arm "RBM Core@C" \
  --total 900 --designs_per_task 100 \
  --num_seq 1 --temp 0.1 --binder_chain_id H --submit

# Assessment & follow-up
python manage_rfa.py assess-rfa-all 6M17 --binder_chain_id H --run_label v1
python manage_rfa.py followup 6M17 --total 2000 --topk 2 --designs_per_task 20 --num_seq 16 --binder_chain_id H --submit
```
Inspect `targets/<PDB>/` for generated artifacts and monitor SLURM job IDs in the
log stream.

## 8. Running the web UI

1. Launch the FastAPI server:
   ```bash
   uvicorn webapp.main:app --reload --port 8000
   ```
2. Open <http://localhost:8000/> in your browser.
3. Use the UI to trigger `init-target`/`decide-scope`/`prep-target`, submit design
   runs, sync assessment outputs, and generate PyMOL or export bundles.
4. When operating with a real cluster, keep the SSH control master active so
   rsync and sbatch invocations remain passwordless.

## 9. Remote filesystem expectations

- Local project data lives under `targets/<PDB>/` with `raw/`, `prep/`, and `designs/`
  subdirectories.
- Cluster-side paths mirror the local layout; ensure sufficient quota for
  `target_root` since AF3 outputs are large.
- Assessment SLURM logs are written under `tools/slurm_logs/` remotely.

Following the steps above yields a workstation that can run the complete
InitBinder workflow locally while coordinating heavy computation on the cluster.
