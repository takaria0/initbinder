"""Run export_files.py via the job store."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path
from typing import List

from .config import load_config
from .job_store import JobStatus, JobStore
from .models import ExportRequest


class ExportError(RuntimeError):
    pass


def _normalize_codon_host(host: str) -> str:
    host_lc = host.strip().lower()
    if host_lc in {"ecoli", "e.coli", "e_coli"}:
        return "e_coli"
    return host_lc


def run_export(request: ExportRequest, *, job_store: JobStore, job_id: str) -> None:
    cfg = load_config()
    workspace = cfg.paths.workspace_root or cfg.paths.project_root
    script = workspace / "export_files.py"
    if not script.exists():
        raise ExportError(f"export_files.py not found at {script}")

    rankings_path = Path(request.rankings_path).expanduser()
    if not rankings_path.exists():
        raise ExportError(f"Rankings TSV not found: {rankings_path}")

    args: List[str] = [
        "python",
        str(script),
        "--rankings_tsv",
        str(rankings_path),
        "--top_n",
        str(request.top_n),
        "--codon_host",
        _normalize_codon_host(request.codon_host),
    ]
    if request.prefix_raw:
        args.extend(["--prefix_raw", request.prefix_raw])
    if request.suffix_raw:
        args.extend(["--suffix_raw", request.suffix_raw])
    if request.gc_target is not None:
        args.extend(["--gc_target", str(request.gc_target)])
    if request.gc_window is not None:
        args.extend(["--gc_window", str(request.gc_window)])
    if request.use_dnachisel:
        args.append("--use_dnachisel")

    env = os.environ.copy()
    existing_path = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = os.pathsep.join(filter(None, [str(cfg.paths.project_root), existing_path]))
    env.setdefault("PYTHONUNBUFFERED", "1")
    env.setdefault("INITBINDER_ROOT", str(workspace))

    logs: List[str] = []

    def _log(line: str) -> None:
        logs.append(line)
        job_store.append_log(job_id, line)

    job_store.update(job_id, status=JobStatus.RUNNING, message="Running export_files.py")
    process = subprocess.Popen(
        args,
        cwd=str(workspace),
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )
    assert process.stdout is not None
    for line in process.stdout:
        _log(line.rstrip())
    rc = process.wait()
    if rc != 0:
        raise ExportError(f"export_files.py exited with code {rc}")

    out_dir = None
    for line in logs:
        if "Exports written to:" in line:
            out_dir = line.split("Exports written to:", 1)[1].strip()
    job_store.update(
        job_id,
        status=JobStatus.SUCCESS,
        message="Export complete",
        details={"out_dir": out_dir, "rankings_path": str(rankings_path)},
    )


__all__ = ["run_export", "ExportError"]
