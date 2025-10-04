"""High-level workflow helpers invoked by API endpoints."""

from __future__ import annotations

import asyncio
from concurrent.futures import ThreadPoolExecutor
from typing import Callable

from .config import load_config
from .job_store import JobStatus, JobStore, get_job_store
from .models import DesignRunRequest, ExportRequest, TargetInitRequest
from .pipeline import PipelineError, init_decide_prep
from .designs import run_design_workflow
from .exporter import ExportError, run_export


_executor: ThreadPoolExecutor | None = None


def _get_executor() -> ThreadPoolExecutor:
    global _executor
    if _executor is None:
        cfg = load_config()
        _executor = ThreadPoolExecutor(max_workers=cfg.background_concurrency)
    return _executor


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
            )
        except PipelineError as exc:
            store.update(job.job_id, status=JobStatus.FAILED, message=str(exc))
        except Exception as exc:  # pragma: no cover - defensive
            store.update(job.job_id, status=JobStatus.FAILED, message=str(exc))

    executor = _get_executor()
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

    executor = _get_executor()
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
    "wait_for_job",
]
