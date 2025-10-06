"""Design pipeline orchestration helpers."""

from __future__ import annotations

import json
import math
import os
import re
import shlex
import textwrap
import threading
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

import yaml

from .config import load_config
from .hpc import ClusterClient
from .job_store import JobStatus, JobStore
from .models import DesignRunRequest

_LAUNCH_PATH_RE = re.compile(r"\[ok] Wrote launcher: (?P<path>.*)")
_JOB_ECHO_RE = re.compile(r"^\[JOB] (?P<var>[A-Za-z0-9_]+)=(?P<jid>\d+)")
_SCRIPT_PATH_RE = re.compile(r"sbatch[^|;]*?([^\s]+\.sh)")


def _parse_launcher_path(log_lines: List[str]) -> Optional[Path]:
    for line in reversed(log_lines):
        match = _LAUNCH_PATH_RE.search(line)
        if match:
            return Path(match.group("path")).expanduser().resolve()
    return None


def _sanitize_epitope_name(name: str) -> str:
    s = str(name).strip()
    return s.replace(" ", "_").replace("/", "_").replace("\\", "_")


def _load_epitope_names(workspace: Path, pdb_id: str) -> List[str]:
    prep_dir = workspace / "targets" / pdb_id.upper() / "prep"
    metadata_path = prep_dir / "epitopes_metadata.json"
    names: List[str] = []
    if metadata_path.exists():
        try:
            data = json.loads(metadata_path.read_text()) or {}
            names = [str(ep.get("name")).strip() for ep in data.get("epitopes", []) if ep.get("name")]
        except Exception:
            names = []
    if not names:
        target_yaml = workspace / "targets" / pdb_id.upper() / "target.yaml"
        if target_yaml.exists():
            try:
                cfg = yaml.safe_load(target_yaml.read_text()) or {}
                names = [str(ep.get("name")).strip() for ep in (cfg.get("epitopes") or []) if ep.get("name")]
            except Exception:
                names = []
    return [name for name in names if name]


def _detect_hotspot_variants(prep_dir: Path, epitope_name: str) -> List[str]:
    san = _sanitize_epitope_name(epitope_name)
    variants: set[str] = set()
    pattern = re.compile(rf"^epitope_{re.escape(san)}_hotspots([A-Za-z0-9]+)\.json$")
    for path in prep_dir.glob(f"epitope_{san}_hotspots*.json"):
        match = pattern.match(path.name)
        if match:
            suffix = match.group(1).upper()
            if suffix:
                variants.add(suffix)
    if not variants:
        if (prep_dir / f"epitope_{san}_hotspots.json").exists():
            variants.add("A")
    if not variants:
        variants.add("A")
    return sorted(variants)


def _discover_arms(workspace: Path, pdb_id: str) -> List[str]:
    prep_dir = workspace / "targets" / pdb_id.upper() / "prep"
    if not prep_dir.exists():
        raise FileNotFoundError(f"Prep directory not found for {pdb_id.upper()} (expected {prep_dir})")

    epitope_names = _load_epitope_names(workspace, pdb_id)
    if not epitope_names:
        raise ValueError("No epitopes available; run prep-target before submitting designs.")

    arms: List[str] = []
    for name in epitope_names:
        variants = _detect_hotspot_variants(prep_dir, name)
        for variant in variants:
            arms.append(f"{name}@{variant}")
    return list(dict.fromkeys(arms))


def _build_pipeline_args(
    request: DesignRunRequest,
    arms: List[str],
    designs_per_task: int,
    run_label: Optional[str],
    binder_chain: str,
) -> List[str]:
    args: List[str] = [request.pdb_id]
    for arm in arms:
        args.extend(["--arm", arm])
    args.extend(["--total", str(request.total_designs)])
    args.extend(["--designs_per_task", str(designs_per_task)])
    args.extend(["--num_seq", str(request.num_sequences)])
    args.extend(["--temp", str(request.temperature)])
    args.extend(["--binder_chain_id", binder_chain])
    if run_label:
        args.extend(["--run_tag", run_label])
    if request.af3_seed is not None:
        args.extend(["--model_seeds", str(request.af3_seed)])
    if request.run_assess:
        args.append("--run_assess")
    return args


def _start_squeue_monitor(
    cluster: ClusterClient,
    job_store: JobStore,
    job_id: str,
    tracked_ids: List[str],
    *,
    user: Optional[str] = None,
    interval_seconds: int = 60,
    max_minutes: int = 360,
) -> None:
    """Spawn a background thread to log `squeue` snapshots for tracked jobs."""

    if not tracked_ids:
        return
    if cluster.cfg.mock:
        job_store.append_log(job_id, "[squeue] Mock cluster; skipping queue monitor.")
        return

    queue_user = (
        user
        or cluster.cfg.user
        or os.getenv("INITBINDER_CLUSTER_USER")
        or os.getenv("USER")
        or "inagakit"
    )

    tracked = [jid for jid in dict.fromkeys(tracked_ids) if jid]
    if not tracked:
        return

    job_store.append_log(
        job_id,
        f"[squeue] Monitoring {len(tracked)} job(s) as {queue_user}: {', '.join(tracked)}",
    )

    interval_seconds = max(15, int(interval_seconds))
    max_minutes = max(5, int(max_minutes))

    def _poll() -> None:
        last_entry: Optional[str] = None
        missing_cycles = 0
        start_ts = time.time()
        while True:
            try:
                job_store.get(job_id)
            except KeyError:
                break

            try:
                result = cluster.squeue(queue_user)
            except Exception as exc:  # pragma: no cover
                job_store.append_log(job_id, f"[squeue][error] {exc}")
                break

            snapshot = (result.stdout or "").strip()
            display = snapshot if snapshot else "(queue empty)"
            stamp = datetime.utcnow().strftime("%H:%M:%S")
            entry = f"[squeue {stamp}]\n{display}"
            if entry != last_entry:
                job_store.append_log(job_id, entry)
                last_entry = entry

            present = any(jid and jid in snapshot for jid in tracked)
            missing_cycles = 0 if present else missing_cycles + 1
            if missing_cycles >= 2:
                job_store.append_log(job_id, "[squeue] Tracked jobs absent for two checks; stopping monitor.")
                break

            elapsed_min = (time.time() - start_ts) / 60.0
            if elapsed_min >= max_minutes:
                job_store.append_log(job_id, f"[squeue] Monitor reached {max_minutes} min limit; stopping.")
                break

            time.sleep(interval_seconds)

    threading.Thread(target=_poll, name=f"squeue-{job_id}", daemon=True).start()


def _schedule_assessment_rescue(
    cluster: ClusterClient,
    job_store: JobStore,
    job_id: str,
    *,
    stage2_ids: List[str],
    pdb_id: str,
    binder_chain: str,
    run_label: str,
) -> None:
    """Submit an assessment job without dependencies if stage2 jobs are stuck.

    The design pipeline launches several stage-2 jobs (AF3 + ProteinMPNN)
    whose successful completion is required before the primary assessment job
    will start.  When any of those stage-2 jobs fail the scheduler leaves the
    downstream assessment job in the ``PD`` (pending) state with a dependency
    related reason (``Dependency`` or ``DependencyNeverSatisfied``).  In that
    situation users typically want to continue with assessment despite the
    partial failure.  This helper polls the scheduler for the stage-2 jobs
    associated with the run.  After three consecutive polls where every job
    remains pending *only* because of dependency issues, it submits a fresh
    assessment job without dependencies.

    The new assessment submission is explicitly marked with
    ``allow_empty_dependencies=True`` so that ``ClusterClient.submit_assessment``
    permits running it even though the dependency list is empty.  All log
    messages are appended to the ``job_store`` so the UI can surface the exact
    decision path to the user.
    """

    if cluster.cfg.mock:
        return
    tracked = [jid for jid in dict.fromkeys(stage2_ids) if jid]
    if not tracked:
        return

    allowed_reasons = {"Dependency", "DependencyNeverSatisfied"}

    def _poll() -> None:
        consecutive = 0
        iterations = 0
        max_checks = 30  # ~1 hour if using default interval
        interval = 120

        while True:
            if iterations >= max_checks:
                job_store.append_log(
                    job_id,
                    "[assessment_rescue] Stage2 dependency stall monitoring timed out.",
                )
                return

            time.sleep(interval)
            iterations += 1

            try:
                snapshot = cluster.describe_jobs(tracked)
            except Exception as exc:  # pragma: no cover - defensive
                job_store.append_log(job_id, f"[assessment_rescue][error] {exc}")
                return

            if not snapshot:
                return

            stalled = all(
                info.get("reason") in allowed_reasons and info.get("state") == "PD"
                for info in snapshot.values()
            )
            if stalled:
                consecutive += 1
            else:
                consecutive = 0
                continue

            if consecutive < 3:
                continue

            job_store.append_log(
                job_id,
                (
                    "[assessment_rescue] Stage2 jobs stalled with dependency issues; "
                    "submitting assessment without dependencies."
                ),
            )

            try:
                rescue_id = cluster.submit_assessment(
                    pdb_id=pdb_id,
                    binder_chain=binder_chain,
                    run_label=run_label,
                    dependencies=[],
                    include_keyword=run_label,
                    allow_empty_dependencies=True,
                )
            except Exception as exc:  # pragma: no cover - defensive
                job_store.append_log(job_id, f"[assessment_rescue][error] {exc}")
                return

            if rescue_id:
                job_store.append_log(
                    job_id,
                    f"[assessment_rescue] Submitted assessment sbatch job {rescue_id}",
                )
                job_store.update(
                    job_id,
                    details={"assessment_rescue_job_id": rescue_id},
                )
                _start_squeue_monitor(cluster, job_store, job_id, [rescue_id])
            return

    threading.Thread(
        target=_poll,
        name=f"assess-rescue-{job_id}",
        daemon=True,
    ).start()


def run_design_workflow(request: DesignRunRequest, *, job_store: JobStore, job_id: str) -> None:
    job_store.update(job_id, status=JobStatus.RUNNING, message="Preparing pipeline configuration")

    cfg = load_config()
    workspace = cfg.paths.workspace_root or cfg.paths.project_root
    arms = _discover_arms(workspace, request.pdb_id)
    if not arms:
        raise ValueError("Failed to determine design arms; ensure prep-target has completed successfully.")

    run_label = request.run_label or datetime.utcnow().strftime("%Y%m%d_%H%M%S")
    designs_per_task = max(1, math.ceil(request.total_designs / max(len(arms), 1)))
    binder_chain = (request.binder_chain_id or "H").strip().upper() or "H"
    assessment_shards = max(1, math.ceil(request.total_designs / 1000))

    job_store.update(
        job_id,
        message=f"Syncing target to cluster ({len(arms)} arms)",
        details={
            "arms": arms,
            "designs_per_task": designs_per_task,
            "run_label": run_label,
            "af3_seed": request.af3_seed,
            "binder_chain": binder_chain,
            "run_assess": request.run_assess,
            "assessment_shards": assessment_shards,
        },
    )
    job_store.append_log(job_id, f"[arms] {', '.join(arms)}")
    job_store.append_log(job_id, f"[designs_per_task] {designs_per_task}")
    job_store.append_log(job_id, f"[af3_seed] {request.af3_seed}")
    job_store.append_log(job_id, f"[binder_chain] {binder_chain}")
    job_store.append_log(job_id, f"[run_assess] {request.run_assess}")
    job_store.append_log(job_id, f"[assessment_shards] {assessment_shards}")

    cluster = ClusterClient(log_hook=lambda line: job_store.append_log(job_id, line))
    local_target_dir = (workspace / "targets" / request.pdb_id.upper()).resolve()
    remote_base = cluster.target_root or cluster.remote_root
    remote_target_desc = (
        str(Path(remote_base) / "targets" / request.pdb_id.upper())
        if remote_base
        else "<remote root unset>"
    )
    job_store.append_log(job_id, f"[rsync] target → {remote_target_desc} (from {local_target_dir})")
    sync_result, backup_rel = cluster.sync_target(request.pdb_id)
    if backup_rel:
        job_store.append_log(job_id, f"[sync] Remote target backup → {backup_rel}")
    job_store.append_log(job_id, f"[cmd] sync_target {request.pdb_id.upper()}")
    if sync_result.stdout:
        for line in sync_result.stdout.splitlines():
            job_store.append_log(job_id, line)
    if sync_result.stderr:
        err = sync_result.stderr.strip()
        if err:
            job_store.append_log(job_id, err)

    # Skip tool sync for now; assume user has done it recently.
    job_store.update(job_id, message="Syncing tools to cluster")
    job_store.append_log(job_id, "[cmd] sync_tools")
    tools_dir = (workspace / "tools").resolve()
    remote_tools_desc = (
        str(Path(cluster.remote_root) / "tools") if cluster.remote_root else "<remote root unset>/tools"
    )
    job_store.append_log(job_id, f"[rsync] tools → {remote_tools_desc} (from {tools_dir})")
    tools_result = cluster.sync_tools()
    if tools_result.stdout:
        for line in tools_result.stdout.splitlines():
            job_store.append_log(job_id, line)
    if tools_result.stderr:
        err = tools_result.stderr.strip()
        if err:
            job_store.append_log(job_id, err)

    pipeline_args = _build_pipeline_args(request, arms, designs_per_task, run_label, binder_chain)
    binder_root = cluster.remote_root
    if binder_root is None:
        raise RuntimeError("Neither target_root nor remote_root configured on cluster; cannot run pipeline")

    env_vars = {"INITBINDER_ROOT": Path(binder_root)}

    target_base: Optional[Path] = None
    if cluster.target_root:
        target_base = Path(cluster.target_root)
    elif cluster.remote_root:
        target_base = Path(cluster.remote_root)

    if target_base is not None:
        # Normalise to the directory that actually contains target folders.
        if target_base.name.lower() != "targets":
            target_base = target_base / "targets"
        env_vars["INITBINDER_TARGET_ROOT"] = target_base

    env_prefix = " ".join(
        f"{key}={shlex.quote(str(value))}" for key, value in env_vars.items()
    )
    remote_cmd = (f"{env_prefix} " if env_prefix else "") + shlex.join([
        "python",
        "manage_rfa.py",
        "pipeline",
        *pipeline_args,
    ])
    job_store.update(job_id, message="Generating pipeline scripts on cluster")
    job_store.append_log(job_id, f"[cmd] {remote_cmd}")
    pipeline_result = cluster.run_in_srun(remote_cmd, use_conda=True)

    log_buffer: List[str] = []
    if pipeline_result.stdout:
        for line in pipeline_result.stdout.splitlines():
            job_store.append_log(job_id, line)
            log_buffer.append(line)
    if pipeline_result.stderr:
        for line in pipeline_result.stderr.splitlines():
            job_store.append_log(job_id, line)
            log_buffer.append(line)

    launcher_path = None
    if log_buffer:
        launcher = _parse_launcher_path(log_buffer)
        if launcher:
            launcher_path = str(launcher)
    if not launcher_path:
        raise FileNotFoundError("Could not locate generated launcher script on the cluster; check manage_rfa output")

    cat_cmd = f"cat {shlex.quote(launcher_path)}"
    launcher_contents = cluster.run(cat_cmd).stdout or ""
    script_order: List[str] = []
    for line in launcher_contents.splitlines():
        match = _SCRIPT_PATH_RE.search(line)
        if match:
            script_order.append(match.group(1))

    job_store.update(job_id, message="Submitting SLURM pipeline")
    submit_cmd = f"conda deactivate >/dev/null 2>&1 || true; conda deactivate >/dev/null 2>&1 || true; bash {shlex.quote(launcher_path)}"
    job_store.append_log(job_id, f"[cmd] {submit_cmd}")
    sbatch_result = cluster.run_in_srun(submit_cmd)
    job_ids: Dict[str, str] = {}
    stage2_ids: List[str] = []
    pipeline_assess_job: Optional[str] = None
    pipeline_assess_merge_job: Optional[str] = None
    if sbatch_result.stdout:
        for line in sbatch_result.stdout.splitlines():
            job_store.append_log(job_id, line)
        idx = 0
        for line in sbatch_result.stdout.splitlines():
            match = re.search(r"Submitted batch job\s+(\d+)", line)
            if match and idx < len(script_order):
                script_path = script_order[idx]
                script_name = Path(script_path).name
                jid = match.group(1)
                job_ids[script_name] = jid
                if "af3batch" in script_name:
                    stage2_ids.append(jid)
                if "assess" in script_name:
                    if "merge" in script_name:
                        pipeline_assess_merge_job = jid
                    else:
                        pipeline_assess_job = jid
                job_store.append_log(job_id, f"[JOB] {script_name}={jid}")
                idx += 1
    if sbatch_result.stderr:
        err = sbatch_result.stderr.strip()
        if err:
            job_store.append_log(job_id, err)

    assessment_job_id: Optional[str] = pipeline_assess_job
    assessment_merge_job_id: Optional[str] = pipeline_assess_merge_job
    needs_manual_assess = bool(stage2_ids) and (not request.run_assess or assessment_job_id is None)
    if needs_manual_assess:
        try:
            fallback_job_id = cluster.submit_assessment(
                pdb_id=request.pdb_id,
                binder_chain=request.binder_chain_id or "H",
                run_label=run_label,
                dependencies=stage2_ids,
                include_keyword=run_label,
            )
        except Exception as exc:  # pragma: no cover
            job_store.append_log(job_id, f"[assessment][error] {exc}")
        else:
            if fallback_job_id:
                assessment_job_id = fallback_job_id
                job_store.append_log(job_id, f"[assessment] Scheduled sbatch job {fallback_job_id}")

    job_store.update(
        job_id,
        details={
            "job_ids": job_ids,
            "remote_launch": launcher_path,
            "assessment_dependencies": stage2_ids,
            "assessment_job_id": assessment_job_id,
            "assessment_via_pipeline": bool(pipeline_assess_job),
            "assessment_merge_job_id": assessment_merge_job_id,
            "run_assess_requested": request.run_assess,
            "assessment_run_label": run_label,
        },
    )

    tracked_for_monitor: List[str] = list(stage2_ids)
    if assessment_job_id:
        tracked_for_monitor.append(assessment_job_id)
    if assessment_merge_job_id:
        tracked_for_monitor.append(assessment_merge_job_id)
    _start_squeue_monitor(cluster, job_store, job_id, tracked_for_monitor)

    if request.run_assess and pipeline_assess_job and stage2_ids:
        # If all stage-2 jobs stall on dependency reasons the rescue monitor
        # will submit a fresh assessment job with no dependencies so users can
        # still review partial results.
        _schedule_assessment_rescue(
            cluster,
            job_store,
            job_id,
            stage2_ids=stage2_ids,
            pdb_id=request.pdb_id,
            binder_chain=binder_chain,
            run_label=run_label,
        )

    job_store.update(
        job_id,
        status=JobStatus.SUCCESS,
        message="Pipeline submitted to cluster",
    )


__all__ = ["run_design_workflow"]
