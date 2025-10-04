"""Local pipeline executors wrapping the existing CLI tools."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path
from typing import Callable, Iterable, List, Optional

from .config import load_config
from .job_store import JobStatus, JobStore

LogCallback = Callable[[str], None]


class PipelineError(RuntimeError):
    pass


def _run_command(cmd: Iterable[str], *, cwd: Path, env: Optional[dict] = None,
                 log: Optional[LogCallback] = None) -> None:
    process_env = os.environ.copy()
    if env:
        process_env.update(env)

    process = subprocess.Popen(
        list(cmd),
        cwd=str(cwd),
        env=process_env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )

    assert process.stdout is not None
    for line in process.stdout:
        if log:
            log(line.rstrip())
    rc = process.wait()
    if rc != 0:
        raise PipelineError(f"Command {' '.join(cmd)} exited with status {rc}")


def run_manage_rfa(subcommand: str, args: List[str], *, log: Optional[LogCallback] = None,
                   env: Optional[dict] = None) -> None:
    cfg = load_config()
    cwd = cfg.paths.workspace_root or cfg.paths.project_root
    script = cwd / "manage_rfa.py"
    if not script.exists():
        raise FileNotFoundError(f"manage_rfa.py not found at {script}")

    cmd = ["python", str(script), subcommand]
    cmd.extend(args)

    # Ensure PYTHONPATH includes project root so modules resolve
    env = env.copy() if env else {}
    existing_path = env.get("PYTHONPATH") or os.environ.get("PYTHONPATH", "")
    env["PYTHONPATH"] = os.pathsep.join(filter(None, [str(cfg.paths.project_root), existing_path]))
    env.setdefault("PYTHONUNBUFFERED", "1")
    env.setdefault("INITBINDER_ROOT", str(cwd))

    _run_command(cmd, cwd=cwd, env=env, log=log)


def init_decide_prep(pdb_id: str, antigen_url: Optional[str], *, job_store: JobStore,
                     job_id: str, run_decide: bool = True, run_prep: bool = True,
                     force: bool = False) -> None:
    def _log(line: str) -> None:
        job_store.append_log(job_id, line)

    job_store.update(job_id, status=JobStatus.RUNNING, message="Initializing target")
    args = [pdb_id]
    if antigen_url:
        args.extend(["--antigen_url", antigen_url])
    if force:
        args.append("--force")
    run_manage_rfa("init-target", args, log=_log)

    if run_decide:
        job_store.update(job_id, message="Running decide-scope")
        run_manage_rfa("decide-scope", [pdb_id], log=_log)
    if run_prep:
        job_store.update(job_id, message="Running prep-target")
        run_manage_rfa("prep-target", [pdb_id], log=_log)

    job_store.update(job_id, status=JobStatus.SUCCESS, message="Target ready")


__all__ = ["PipelineError", "run_manage_rfa", "init_decide_prep"]
