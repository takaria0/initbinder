# Design Submission Flow Notes

This document walks through what happens after a user presses **Submit to cluster** in the InitBinder UI. It captures the current codepaths for both the **RFantibody** and **BoltzGen** design engines, the intermediate data structures that are produced, and the commands ultimately executed on the remote cluster.

---

## 1. Frontend → API

1. The design form lives in `webapp/templates/index.html`. Engine-specific labels and visibility are managed in `webapp/static/js/app.js`.
2. Clicking **Submit to cluster** triggers `queueDesignRun()` (app.js:2135).
   - The handler builds a JSON payload containing:
     - `pdb_id` – current target (`state.currentPdb`).
     - `model_engine` – selected engine (`state.selectedDesignEngine`).
     - Shared fields that remain visible for the active engine (e.g. `total_designs`, `run_label`, `run_assess`, `rfdiff_crop_radius`).
     - RFantibody-only fields (`num_sequences`, `temperature`, `binder_chain_id`, `af3_seed`) are included only when that engine marks them as “active”.
   - The payload is posted to `POST /api/designs/run`.
3. The UI disables the submit button, resets the job log panel, and starts polling once a job ID is returned.

Example RFantibody payload:

```json
{
  "pdb_id": "6M17",
  "model_engine": "rfantibody",
  "total_designs": 900,
  "num_sequences": 1,
  "temperature": 0.1,
  "binder_chain_id": "H",
  "af3_seed": 1,
  "run_label": "demo_run",
  "run_assess": true,
  "rfdiff_crop_radius": null
}
```

Example BoltzGen payload (note the missing RFantibody-only fields):

```json
{
  "pdb_id": "6M17",
  "model_engine": "boltzgen",
  "total_designs": 1000,
  "run_label": "bg_test",
  "run_assess": true
}
```

---

## 2. FastAPI Layer

1. `webapp/main.py` registers `api_design_run()` and returns the `job_id` created by the workflow layer.
2. The request object is validated by `DesignRunRequest` (Pydantic model in `webapp/models.py`) which ensures correct types and default values.

---

## 3. Workflow Submission (`webapp/workflows.py`)

1. `submit_design_run()` is called with the validated `DesignRunRequest`.
   - A new job entry is created in the job store (`design_run` kind). The label uses the upper-cased model engine (`RFANTIBODY` or `BOLTZGEN`) and the PDB ID.
   - The request payload is stored in the job details for later debugging.
   - Work is scheduled on the shared background executor (`ThreadPoolExecutor`).
2. The worker thread calls `run_design_workflow()` in `webapp/designs.py` with the job context.

---

## 4. Design Engine Registry (`webapp/designs.py`)

1. `run_design_workflow()` looks up the configured engine in the registry:
   - `RFAntibodyEngine` (default, RFdiffusion → ProteinMPNN → AlphaFold3 workflow).
   - `BoltzGenEngine` (BoltzGen diffusion pipeline).
2. Each engine implementation:
   - Logs `model_engine` to the job log.
   - Updates the job message/status as it moves through the pipeline.
   - Emits detailed log entries so the UI can surface progress.

### 4.1 RFantibody Engine Flow

1. **Arm Discovery**
   - `_discover_arms()` collects epitope variants from the local workspace (`targets/<PDB>/prep/`).
   - Arms (e.g. `"Receptor Binding Motif Core@A"`) are logged and stored in the job details.
2. **Cluster Sync**
   - `ClusterClient.sync_target(pdb_id)` copies the target folder to the remote cluster, including epitope metadata and prepared structures.
   - `ClusterClient.sync_tools()` syncs the `tools/` directory (includes launcher templates).
3. **Pipeline Script Generation**
   - Builds the CLI arguments for `manage_rfa.py pipeline` (e.g., `--arm`, `--total`, `--designs_per_task`, `--num_seq`, `--temp`, `--model_seeds`, `--binder_chain_id`, `--run_tag`).
   - Executes on the cluster via `cluster.run_in_srun(remote_cmd, use_conda=True)`.
   - The command has the general shape:
     ```
     INITBINDER_ROOT=<remote_root> INITBINDER_TARGET_ROOT=<remote_targets>
     python manage_rfa.py pipeline <PDB> --arm ... --total ... --designs_per_task ... --num_seq ...
     ```
   - The `manage_rfa.py` execution writes sbatch scripts and a launcher shell script on the cluster; its log is streamed back and stored.
4. **Launcher Extraction**
   - Parses the `manage_rfa` logs to find the path of the generated launcher (`tools/launchers/launch_pipeline_<pdb>_<timestamp>.sh`).
   - Reads the launcher to capture the order of sbatch submissions (RFdiffusion arrays, ProteinMPNN jobs, AlphaFold3 stage 1/2, assessment shards).
5. **Pipeline Submission**
   - Runs the launcher via `cluster.run_in_srun("bash <launcher>")`.
   - Captures the SLURM job IDs from stdout and records them in the job store (`assessment_dependencies`, `assessment_job_id`, etc.).
   - Starts an `squeue` monitor thread to stream queue status into the job log.
6. **Assessment Handling**
   - If the user enabled `run_assess`, monitors for stuck stage-2 jobs and can submit an assessment rescue job without dependencies.
   - Starts a background thread to poll for AF3 rankings on the cluster and triggers `sync_assessments_back` when ready.

Resulting remote commands (RFantibody):

```
INITBINDER_ROOT=/data/... INITBINDER_TARGET_ROOT=/pub/.../targets \
python manage_rfa.py pipeline 6M17 \
  --arm "Receptor Binding Motif Core@A" \
  --arm "Receptor Binding Motif Core@B" \
  --total 900 \
  --designs_per_task 100 \
  --num_seq 1 \
  --temp 0.1 \
  --binder_chain_id H \
  --model_seeds 1 \
  --run_tag demo_run

bash tools/launchers/launch_pipeline_6M17_20251011_221530.sh
```

### 4.2 BoltzGen Engine Flow

1. **Spec Generation**
   - Creates a workspace subdirectory: `targets/<PDB>/designs/_boltzgen/<run_label>/`.
   - Writes a single BoltzGen design YAML that references the prepared PDB without epitope/binding annotations. Any detected prep arms are logged for context but ignored so the requested total designs flow directly into the run.
   - Logs spec paths and hotspot counts, updating job details (`spec_paths`, `spec_outputs`).
2. **Cluster Sync**
   - Same target and tools rsync steps as RFantibody (ensures remote has the fresh prep data and tools/ scripts).
3. **BoltzGen Pipeline Command**
   - Builds a remote command for the helper script `tools/boltzgen/pipeline.py`.
   - Injects resource configuration from `cfg/webapp.yaml` (`cluster.boltzgen`) such as partition, account, GPUs, CPUs, memory, walltime, cache directory, and conda activation string.
   - The command looks like:
     ```
     INITBINDER_ROOT=/data/... INITBINDER_TARGET_ROOT=/pub/.../targets \
     python tools/boltzgen/pipeline.py pipeline 6M17 \
       --run_label bg_test \
       --num_designs 1000 \
       --protocol protein-anything \
       --scripts_dir tools/boltzgen \
       --launcher_dir tools/launchers \
       --output_root /pub/.../targets/6M17/designs/_boltzgen/bg_test \
       --partition gpu \
       --account ccl_lab_gpu \
       --gpus A100:1 \
       --cpus 8 \
     --mem 64G \
     --time_h 12 \
     --conda_activate "conda activate boltzgen" \
     --spec /pub/.../targets/6M17/designs/_boltzgen/bg_test/full_target.yaml \
     --submit
     ```
4. **Launcher Output**
   - `tools/boltzgen/pipeline.py` writes per-spec sbatch scripts under `tools/boltzgen/` and a launcher under `tools/launchers/`.
   - Submits the launcher immediately (`--submit` behaviour), capturing job IDs from stdout.
   - Parses the helper output to record the launcher path, submitted job IDs, and the final run label.
5. **Monitoring**
   - Uses the same `squeue` monitoring as RFantibody, logging cluster status updates.
   - Currently `run_assess` is logged as not yet supported (assessment jobs are skipped until a downstream integration exists).

Resulting remote commands (BoltzGen):

```
INITBINDER_ROOT=/data/... INITBINDER_TARGET_ROOT=/pub/.../targets \
python tools/boltzgen/pipeline.py pipeline 6M17 \
  --run_label bg_test \
  --num_designs 1000 \
  --protocol protein-anything \
  --scripts_dir tools/boltzgen \
  --launcher_dir tools/launchers \
  --output_root /pub/.../targets/6M17/designs/_boltzgen/bg_test \
  --partition gpu \
  --account ccl_lab_gpu \
  --gpus A100:1 \
  --cpus 8 \
  --mem 64G \
  --time_h 12 \
  --conda_activate "conda activate boltzgen" \
  --spec /pub/.../targets/6M17/designs/_boltzgen/bg_test/full_target.yaml \
  --submit
```

---

## 5. Cluster Helper Script (`tools/boltzgen/pipeline.py`)

1. Builds sbatch scripts for each spec with the following structure:
   - SLURM headers populated from CLI arguments (`--partition`, `--account`, `--gpus`, `--cpus`, `--mem`, `--time_h`).
   - Optional conda activation snippet for the BoltzGen environment.
   - `boltzgen run <spec> --output <dir> --protocol ... --num_designs ...`.
   - Optional Hugging Face cache override if `--cache_dir` is provided.
2. Generates a launcher that simply `sbatch`s all scripts sequentially.
3. With `--submit`, the launcher is executed immediately; submitted job IDs are printed to stdout for upstream logging.

## 5.1 BoltzGen Environment on the Cluster

- The helper script activates whatever environment is referenced by `cluster.boltzgen.conda_activate` in `cfg/webapp.yaml` (default: `conda activate boltzgen`).
- Recommended location: create a dedicated environment under your cluster home, e.g. `~/miniconda3/envs/boltzgen` (resolves to `/data/homezvol1/inagakit/miniconda3/envs/boltzgen` on HPC3).
- Installation steps on the cluster login node:
  ```bash
  module load conda  # if required by your site setup
  conda create -n boltzgen python=3.11
  conda activate boltzgen
  pip install boltzgen
  ```
- Ensure the environment provides the `boltzgen` CLI and GPU-aware dependencies (PyTorch, etc.); the pip package pulls in all runtime requirements.
- If you relocate the environment or use a module system, update `cluster.boltzgen.conda_activate` accordingly (e.g. `source /opt/boltzgen/env.sh`).
- Large model downloads default to `~/.cache` on the cluster. To steer them elsewhere (e.g. fast scratch), set `cluster.boltzgen.cache_dir` in `cfg/webapp.yaml`; the helper automatically exports `HF_HOME=<cache_dir>`.

---

## 6. Viewing BoltzGen Results in the UI

1. **Sync outputs** – Click *Sync BoltzGen results* (or call `POST /api/targets/{PDB}/boltzgen/sync`) to rsync `targets/<PDB>/designs/_boltzgen/<run_label>` back from the cluster. Pass `run_label` to limit the transfer, otherwise the entire `_boltzgen/` tree is mirrored.
2. **Choose the source** – In the Results panel select *BoltzGen metrics*. A dedicated run list appears below the AF3 history, showing every run/spec combo along with a ✔ if `final_ranked_designs/all_designs_metrics.csv` exists locally.
3. **Load metrics** – Enter the run label (and optional spec name) and press *Refresh*. The backend parses the `all_designs_metrics.csv` file, maps `design_to_target_iptm → ipTM` and `filter_rmsd → RMSD`, and feeds the rows into the existing table + scatter plot so you can inspect ipTM vs RMSD immediately.
4. **Feature scope** – AF3-specific tools (plot generation, Golden Gate planner, PyMOL launchers, exports) stay disabled while the BoltzGen source is active. The sortable table, scatter plot, metadata drawer, and job history continue to function.

| Handy reference | Path / Action |
| --- | --- |
| Cluster metrics CSV | `/pub/.../targets/<PDB>/designs/_boltzgen/<run>/<spec>/final_ranked_designs/all_designs_metrics.csv` |
| Local mirror | `<repo>/targets/<PDB>/designs/_boltzgen/<run>/<spec>/final_ranked_designs/all_designs_metrics.csv` |
| UI controls | Results → Source dropdown + *Sync BoltzGen results* |

## 7. Job Store & UI Updates

1. All engines write log lines via `job_store.append_log(job_id, message)` (FastAPI streams them to the UI via polling).
2. Additional metadata is stored in job details for downstream panels:
   - `arms`, `designs_per_task`, `spec_paths`, `job_ids`, `remote_launch`, `boltzgen_output_root`, `assessment_*`.
3. The UI polls `/api/jobs/{job_id}` to display progress, render job logs, and enable follow-up actions (syncing assessments, opening PyMOL bundles, etc.).

---

## 8. Quick Reference of Key Modules

| Component | Responsibility |
| --- | --- |
| `webapp/static/js/app.js` | UI behaviour, engine selector, request payload assembly |
| `webapp/models.py` | Pydantic schemas for design requests/responses, engine metadata |
| `webapp/main.py` | FastAPI endpoints bridging UI to workflows |
| `webapp/workflows.py` | Job creation, executor management, invocation of engines |
| `webapp/designs.py` | Engine registry, RFantibody + BoltzGen implementations |
| `webapp/hpc.py` | ClusterClient: rsync, ssh, sbatch, squeue helpers |
| `manage_rfa.py` | Generates RFdiffusion/ProteinMPNN/AF3 sbatch scripts on cluster |
| `tools/boltzgen/pipeline.py` | Generates BoltzGen sbatch scripts/launcher and submits them |

---

## 9. Summary

- The GUI simply posts structured JSON to a single endpoint. All pipeline specificity is handled server-side via pluggable engines.
- RFantibody continues to drive the legacy `manage_rfa.py` pipeline, while BoltzGen uses new YAML specs and a concise cluster helper script.
- The job store maintains full transparency: generated files, launchers, job IDs, and command lines are logged for debugging and replay.
- Future engines can be integrated by implementing a new `DesignEngine` class, defining appropriate UI field metadata, and registering it with the engine registry. The frontend will automatically pick up labels, visibility, and descriptions from the backend payload. 
