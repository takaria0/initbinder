"""Design pipeline orchestration helpers."""

from __future__ import annotations

import re
import shlex
from pathlib import Path
from typing import Dict, List, Optional

from .config import load_config
from .hpc import ClusterClient
from .job_store import JobStatus, JobStore
from .models import DesignRunRequest
from .pipeline import run_manage_rfa

_LAUNCH_PATH_RE = re.compile(r"\[ok] Wrote launcher: (?P<path>.*)")
_JOB_ECHO_RE = re.compile(r"^\[JOB] (?P<var>[A-Za-z0-9_]+)=(?P<jid>\d+)")


def _instrument_launcher(script_path: Path) -> None:
    lines = script_path.read_text().splitlines()
    patched: List[str] = []
    assign_pattern = re.compile(r'^([A-Za-z0-9_]+)=\$\(\s*sbatch\b')
    for line in lines:
        patched.append(line)
        stripped = line.strip()
        m = assign_pattern.match(stripped)
        if m:
            var_name = m.group(1)
            patched.append(f'echo "[JOB] {var_name}=${{{var_name}}}"')
    script_path.write_text("\n".join(patched) + "\n")


def _parse_launcher_path(log_lines: List[str]) -> Optional[Path]:
    for line in reversed(log_lines):
        match = _LAUNCH_PATH_RE.search(line)
        if match:
            return Path(match.group("path")).expanduser().resolve()
    return None


def _build_pipeline_args(request: DesignRunRequest) -> List[str]:
    args: List[str] = [request.pdb_id]
    for arm in request.arms:
        arm_clean = arm.strip()
        if arm_clean:
            args.extend(["--arm", arm_clean])
    if request.total_designs:
        args.extend(["--total", str(request.total_designs)])
    if request.designs_per_task:
        args.extend(["--designs_per_task", str(request.designs_per_task)])
    if request.num_sequences:
        args.extend(["--num_seq", str(request.num_sequences)])
    if request.temperature is not None:
        args.extend(["--temp", str(request.temperature)])
    if request.binder_chain_id:
        args.extend(["--binder_chain_id", request.binder_chain_id])
    if request.run_label:
        args.extend(["--run_tag", request.run_label])
    return args


def run_design_workflow(request: DesignRunRequest, *, job_store: JobStore, job_id: str) -> None:
    job_store.update(job_id, status=JobStatus.RUNNING, message="Generating pipeline scripts")
    log_buffer: List[str] = []

    def _log(line: str) -> None:
        log_buffer.append(line)
        job_store.append_log(job_id, line)

    args = _build_pipeline_args(request)
    run_manage_rfa("pipeline", args, log=_log)

    launcher = _parse_launcher_path(log_buffer)
    if not launcher or not launcher.exists():
        raise FileNotFoundError("Could not locate generated launcher script; check manage_rfa output")

    _instrument_launcher(launcher)
    job_store.update(job_id, message=f"Launcher instrumented: {launcher}")

    cfg = load_config()
    workspace = cfg.paths.workspace_root or cfg.paths.project_root
    try:
        rel_path = launcher.relative_to(workspace)
    except ValueError as exc:
        raise RuntimeError(f"Launcher {launcher} is outside workspace root {workspace}") from exc

    cluster = ClusterClient()
    job_store.update(job_id, message="Syncing target folder to cluster")
    cluster.sync_target(request.pdb_id)
    job_store.append_log(job_id, "[sync] target directory copied")
    job_store.update(job_id, message="Syncing tools to cluster")
    cluster.sync_tools()
    job_store.append_log(job_id, "[sync] tools directory copied")

    remote_path = cluster.remote_path(rel_path)
    job_store.update(job_id, message=f"Submitting SLURM pipeline via {remote_path}")
    remote_cmd = f"bash {shlex.quote(str(remote_path))}"
    result = cluster.ssh(remote_cmd)
    if result.stderr:
        job_store.append_log(job_id, result.stderr)
    if result.stdout:
        for line in result.stdout.splitlines():
            job_store.append_log(job_id, line)
    job_ids: Dict[str, str] = {}
    for line in result.stdout.splitlines():
        match = _JOB_ECHO_RE.match(line.strip())
        if match:
            job_ids[match.group("var")] = match.group("jid")
    job_store.update(
        job_id,
        status=JobStatus.SUCCESS,
        message="Pipeline submitted to cluster",
        details={"job_ids": job_ids, "remote_launch": str(remote_path)},
    )


__all__ = ["run_design_workflow"]
