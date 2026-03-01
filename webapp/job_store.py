"""In-memory + optional on-disk job status tracking for long-running backend tasks."""

from __future__ import annotations

import json
import threading
import time
import uuid
from dataclasses import asdict, dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional


class JobStatus(str, Enum):
    PENDING = "pending"
    RUNNING = "running"
    SUCCESS = "success"
    FAILED = "failed"
    CANCELED = "canceled"


@dataclass
class JobRecord:
    job_id: str
    kind: str
    label: str
    status: JobStatus = JobStatus.PENDING
    created_at: float = field(default_factory=time.time)
    started_at: Optional[float] = None
    finished_at: Optional[float] = None
    progress: Optional[float] = None
    message: Optional[str] = None
    details: Dict[str, Any] = field(default_factory=dict)
    logs: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        payload = asdict(self)
        payload["status"] = self.status.value
        return payload


class JobStore:
    def __init__(self, persist_dir: Optional[Path] = None) -> None:
        self._lock = threading.RLock()
        self._jobs: Dict[str, JobRecord] = {}
        self._persist_dir = persist_dir
        if self._persist_dir:
            self._persist_dir.mkdir(parents=True, exist_ok=True)
            self._restore_from_disk()

    # -- job lifecycle -----------------------------------------------------
    def create_job(self, kind: str, label: str, *, job_id: Optional[str] = None,
                   initial_status: JobStatus = JobStatus.PENDING,
                   details: Optional[Dict[str, Any]] = None) -> JobRecord:
        with self._lock:
            job_id = job_id or uuid.uuid4().hex
            record = JobRecord(job_id=job_id, kind=kind, label=label,
                               status=initial_status,
                               details=details.copy() if details else {})
            self._jobs[job_id] = record
            self._flush(record)
            return record

    def update(self, job_id: str, *, status: Optional[JobStatus] = None,
               message: Optional[str] = None, progress: Optional[float] = None,
               append_log: Optional[str] = None, details: Optional[Dict[str, Any]] = None) -> JobRecord:
        with self._lock:
            record = self._jobs[job_id]
            if status is not None:
                record.status = status
                if status == JobStatus.RUNNING and record.started_at is None:
                    record.started_at = time.time()
                if status in {JobStatus.SUCCESS, JobStatus.FAILED, JobStatus.CANCELED}:
                    record.finished_at = time.time()
            if message is not None:
                record.message = message
            if progress is not None:
                record.progress = max(0.0, min(1.0, progress))
            if append_log:
                record.logs.append(append_log)
            if details:
                record.details.update(details)
            self._flush(record)
            return record

    def append_log(self, job_id: str, line: str) -> None:
        self.update(job_id, append_log=line)

    def get(self, job_id: str) -> JobRecord:
        with self._lock:
            return self._jobs[job_id]

    def list_jobs(self) -> List[JobRecord]:
        with self._lock:
            return list(self._jobs.values())

    def _flush(self, record: JobRecord) -> None:
        if not self._persist_dir:
            return
        path = self._persist_dir / f"{record.job_id}.json"
        path.write_text(json.dumps(record.to_dict(), indent=2))

    # -- persistence helpers -----------------------------------------------
    def _restore_from_disk(self) -> None:
        if not self._persist_dir:
            return

        candidates = sorted(self._persist_dir.glob("*.json"))
        restored: Dict[str, JobRecord] = {}
        for path in candidates:
            try:
                payload = json.loads(path.read_text())
            except Exception:
                continue

            if not isinstance(payload, dict):
                continue

            def _coerce_float(value: Optional[object]) -> Optional[float]:
                if value is None:
                    return None
                try:
                    return float(value)
                except (TypeError, ValueError):
                    return None

            job_id = str(payload.get("job_id") or path.stem)
            kind = str(payload.get("kind") or "unknown")
            label = str(payload.get("label") or job_id)
            status_raw = payload.get("status", JobStatus.PENDING.value)
            try:
                status = JobStatus(status_raw)
            except ValueError:
                status = JobStatus.PENDING

            record = JobRecord(
                job_id=job_id,
                kind=kind,
                label=label,
                status=status,
                created_at=_coerce_float(payload.get("created_at")) or path.stat().st_mtime,
                started_at=_coerce_float(payload.get("started_at")),
                finished_at=_coerce_float(payload.get("finished_at")),
                progress=_coerce_float(payload.get("progress")),
                message=payload.get("message"),
                details=payload.get("details") or {},
                logs=list(payload.get("logs") or []),
            )
            restored[job_id] = record

        if not restored:
            return

        # Limit to most recent entries by created_at to keep memory bounded.
        recent = sorted(restored.values(), key=lambda rec: rec.created_at, reverse=True)
        for record in recent[:500]:
            self._jobs[record.job_id] = record


_default_store: Optional[JobStore] = None


def get_job_store(persist_dir: Optional[Path] = None) -> JobStore:
    global _default_store
    if _default_store is None:
        _default_store = JobStore(persist_dir=persist_dir)
    return _default_store


__all__ = [
    "JobStatus",
    "JobRecord",
    "JobStore",
    "get_job_store",
]
