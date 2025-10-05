"""High-level workflow helpers invoked by API endpoints."""

from __future__ import annotations

import asyncio
import re
import shlex
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Callable

from . import preferences
from .config import load_config
from .hpc import ClusterClient
from .job_store import JobStatus, JobStore, get_job_store
from .models import (
    AssessmentRunRequest,
    DesignRunRequest,
    ExportRequest,
    TargetInitRequest,
)
from .pipeline import PipelineError, init_decide_prep
from .designs import run_design_workflow
from .exporter import ExportError, run_export


_executor: ThreadPoolExecutor | None = None
_init_executor: ThreadPoolExecutor | None = None
_export_executor: ThreadPoolExecutor | None = None
_cluster_executor: ThreadPoolExecutor | None = None


def _get_executor() -> ThreadPoolExecutor:
    global _executor
    if _executor is None:
        cfg = load_config()
        _executor = ThreadPoolExecutor(max_workers=cfg.background_concurrency)
    return _executor


def _get_init_executor() -> ThreadPoolExecutor:
    global _init_executor
    if _init_executor is None:
        cfg = load_config()
        workers = max(1, min(2, cfg.background_concurrency or 1))
        _init_executor = ThreadPoolExecutor(max_workers=workers)
    return _init_executor


def _get_export_executor() -> ThreadPoolExecutor:
    """Return a lightweight executor for short-lived export jobs."""
    global _export_executor
    if _export_executor is None:
        cfg = load_config()
        # Ensure at least one worker even if background_concurrency == 0.
        # Limit to a small number so exports cannot overwhelm the system.
        workers = max(1, min(2, cfg.background_concurrency or 1))
        _export_executor = ThreadPoolExecutor(max_workers=workers)
    return _export_executor


def _get_cluster_executor() -> ThreadPoolExecutor:
    """Dedicated executor for cluster sync operations."""
    global _cluster_executor
    if _cluster_executor is None:
        cfg = load_config()
        workers = max(1, min(2, cfg.background_concurrency or 1))
        _cluster_executor = ThreadPoolExecutor(max_workers=workers)
    return _cluster_executor


def submit_target_initialization(request: TargetInitRequest, *,
                                  job_store: JobStore | None = None) -> str:
    store = job_store or get_job_store(load_config().log_dir)
    label = f"Init target {request.pdb_id.upper()}"
    job = store.create_job("target_init", label, details=request.dict())

    def _run() -> None:
        try:
            init_decide_prep(
                request.pdb_id,
                request.antigen_url,
                job_store=store,
                job_id=job.job_id,
                run_decide=request.run_decide_scope,
                run_prep=request.run_prep,
                force=request.force_refresh,
                num_epitopes=request.num_epitopes,
            )
            try:
                preferences.record_target_usage(
                    pdb_id=request.pdb_id,
                    name=request.preset_name,
                    antigen_url=request.antigen_url,
                    num_epitopes=request.num_epitopes,
                )
            except Exception as exc:  # pragma: no cover - preferences are best effort
                store.append_log(job.job_id, f"[warn] failed to update target presets: {exc}")
        except PipelineError as exc:
            store.update(job.job_id, status=JobStatus.FAILED, message=str(exc))
        except Exception as exc:  # pragma: no cover - defensive
            store.update(job.job_id, status=JobStatus.FAILED, message=str(exc))

    executor = _get_init_executor()
    executor.submit(_run)
    return job.job_id


def submit_design_run(request: DesignRunRequest, *, job_store: JobStore | None = None) -> str:
    store = job_store or get_job_store(load_config().log_dir)
    label = f"Design pipeline {request.pdb_id.upper()}"
    job = store.create_job("design_run", label, details=request.dict())

    def _run() -> None:
        try:
            run_design_workflow(request, job_store=store, job_id=job.job_id)
        except Exception as exc:  # pragma: no cover - defensive
            store.update(job.job_id, status=JobStatus.FAILED, message=str(exc))

    executor = _get_executor()
    executor.submit(_run)
    return job.job_id


def submit_export(request: ExportRequest, *, job_store: JobStore | None = None) -> str:
    store = job_store or get_job_store(load_config().log_dir)
    label = f"Export designs {request.pdb_id.upper()}"
    job = store.create_job("export", label, details=request.dict())

    def _run() -> None:
        try:
            run_export(request, job_store=store, job_id=job.job_id)
        except ExportError as exc:
            store.update(job.job_id, status=JobStatus.FAILED, message=str(exc))
        except Exception as exc:  # pragma: no cover - defensive
            store.update(job.job_id, status=JobStatus.FAILED, message=str(exc))

    executor = _get_export_executor()
    executor.submit(_run)
    return job.job_id


def submit_assessment_run(request: AssessmentRunRequest, *, job_store: JobStore | None = None) -> str:
    store = job_store or get_job_store(load_config().log_dir)
    label = f"Assess designs {request.pdb_id.upper()}"
    job = store.create_job("assessment", label, details=request.dict())

    def _run() -> None:
        include_keyword = request.include_keyword or request.run_label
        args = [
            request.pdb_id,
            "--binder_chain_id",
            request.binder_chain_id,
            "--run_label",
            request.run_label,
            "--sbatch",
            "--submit",
        ]
        if include_keyword:
            args.extend(["--include_keyword", include_keyword])

        client = ClusterClient()
        store.update(job.job_id, status=JobStatus.RUNNING, message="Submitting assessment jobs")

        try:
            store.append_log(job.job_id, "[cmd] sync_tools")
            sync_result = client.sync_tools()
            if sync_result.stdout:
                for line in sync_result.stdout.splitlines():
                    store.append_log(job.job_id, line)
            if sync_result.stderr:
                err = sync_result.stderr.strip()
                if err:
                    store.append_log(job.job_id, err)
        except Exception as exc:
            store.append_log(job.job_id, f"[warn] sync_tools failed: {exc}")

        remote_root = client.remote_root
        if remote_root is None:
            store.update(job.job_id, status=JobStatus.FAILED, message="Cluster remote_root not configured")
            return

        cmd = ["python", "manage_rfa.py", "assess-rfa-all", *args]
        env_prefix = f"INITBINDER_ROOT={shlex.quote(str(remote_root))}"
        remote_cmd = env_prefix + " " + shlex.join(cmd)
        store.append_log(job.job_id, f"[cmd] {remote_cmd}")

        try:
            result = client.run(remote_cmd, use_conda=True)
            job_ids: list[str] = []
            if result.stdout:
                for line in result.stdout.splitlines():
                    store.append_log(job.job_id, line)
                    match = re.search(r"Submitted batch job (\d+)", line)
                    if match:
                        job_ids.append(match.group(1))
            if result.stderr:
                for line in result.stderr.splitlines():
                    store.append_log(job.job_id, line)

            store.update(
                job.job_id,
                status=JobStatus.SUCCESS,
                message="Assessment submitted to cluster",
                details={"sbatch_job_ids": job_ids, "run_label": request.run_label},
            )
        except Exception as exc:  # pragma: no cover - defensive
            store.update(job.job_id, status=JobStatus.FAILED, message=str(exc))

    executor = _get_executor()
    executor.submit(_run)
    return job.job_id


def submit_assessment_sync(pdb_id: str, run_label: str | None = None,
                           *, job_store: JobStore | None = None) -> str:
    store = job_store or get_job_store(load_config().log_dir)
    label = f"Sync assessments {pdb_id.upper()}"
    if run_label:
        label += f" ({run_label})"
    details = {"pdb_id": pdb_id, "run_label": run_label}
    job = store.create_job("assessment_sync", label, details=details)

    def _run() -> None:
        client = ClusterClient()
        store.update(job.job_id, status=JobStatus.RUNNING, message="Syncing assessments from cluster")

        base_rel = Path("targets") / pdb_id.upper() / "designs"
        if run_label:
            rel = base_rel / "_assessments" / run_label
        else:
            rel = base_rel
        local_path = (client.local_root / rel).resolve()
        remote_root = client.target_root or client.remote_root
        remote_path = str(Path(remote_root) / rel) if remote_root else None

        try:
            result = client.sync_assessments_back(pdb_id, run_label=run_label)
            if result.stdout:
                for line in result.stdout.splitlines():
                    store.append_log(job.job_id, line)
            if result.stderr:
                for line in result.stderr.splitlines():
                    store.append_log(job.job_id, line)

            message = "Assessments already synced" if result.skipped else "Assessments synced"
            store.update(
                job.job_id,
                status=JobStatus.SUCCESS,
                message=message,
                details={
                    "pdb_id": pdb_id,
                    "run_label": run_label,
                    "local_path": str(local_path),
                    "remote_path": str(remote_path) if remote_path else None,
                    "skipped": result.skipped,
                },
            )
        except Exception as exc:  # pragma: no cover - defensive
            store.update(job.job_id, status=JobStatus.FAILED, message=str(exc))

    executor = _get_cluster_executor()
    executor.submit(_run)
    return job.job_id


async def wait_for_job(job_id: str, *, job_store: JobStore | None = None,
                       poll_interval: float = 0.5) -> None:
    store = job_store or get_job_store()
    while True:
        record = store.get(job_id)
        if record.status not in {JobStatus.PENDING, JobStatus.RUNNING}:
            break
        await asyncio.sleep(poll_interval)


__all__ = [
    "submit_target_initialization",
    "submit_design_run",
    "submit_export",
    "submit_assessment_run",
    "submit_assessment_sync",
    "wait_for_job",
]
