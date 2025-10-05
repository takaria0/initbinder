"""Local pipeline executors wrapping the existing CLI tools."""

from __future__ import annotations

import os
import shutil
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Callable, Iterable, List, Optional

import yaml

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


def _persist_num_epitopes(pdb_id: str, num_epitopes: int) -> None:
    cfg = load_config()
    targets_dir = cfg.paths.targets_dir or (cfg.paths.workspace_root / "targets")
    target_dir = targets_dir / pdb_id.upper()
    target_yaml = target_dir / "target.yaml"
    if not target_yaml.exists():
        return

    try:
        data = yaml.safe_load(target_yaml.read_text()) or {}
    except Exception:
        data = {}

    webapp_block = data.setdefault("webapp", {})
    webapp_block["num_epitopes"] = int(num_epitopes)

    try:
        target_yaml.write_text(yaml.safe_dump(data, sort_keys=False))
    except Exception as exc:
        raise PipelineError(f"Failed to persist num_epitopes to {target_yaml}: {exc}") from exc


def init_decide_prep(
    pdb_id: str,
    antigen_url: Optional[str],
    *,
    job_store: JobStore,
    job_id: str,
    run_decide: bool = True,
    run_prep: bool = True,
    force: bool = False,
    num_epitopes: Optional[int] = None,
) -> None:
    def _log(line: str) -> None:
        job_store.append_log(job_id, line)

    job_store.update(job_id, status=JobStatus.RUNNING, message="Initializing target")
    snapshot_dir = _snapshot_target_state(pdb_id)
    if snapshot_dir:
        job_store.append_log(job_id, f"[snapshot] Saved prior target state → {snapshot_dir}")
        job_store.update(job_id, details={"snapshot_dir": str(snapshot_dir)})
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

    if num_epitopes is not None:
        try:
            _persist_num_epitopes(pdb_id, num_epitopes)
            job_store.append_log(job_id, f"[prefs] recorded num_epitopes={num_epitopes}")
        except Exception as exc:  # pragma: no cover - best effort to persist metadata
            job_store.append_log(job_id, f"[warn] failed to persist num_epitopes: {exc}")

    job_store.update(job_id, status=JobStatus.SUCCESS, message="Target ready")


def _snapshot_target_state(pdb_id: str, *, keep: int = 8) -> Optional[Path]:
    cfg = load_config()
    targets_dir = cfg.paths.targets_dir or (cfg.paths.workspace_root / "targets")
    target_dir = (targets_dir / pdb_id.upper()).resolve()
    if not target_dir.exists():
        return None

    target_yaml = target_dir / "target.yaml"
    prep_dir = target_dir / "prep"
    if not target_yaml.exists() and not prep_dir.exists():
        return None

    history_root = target_dir / "history"
    timestamp = datetime.utcnow().strftime("%Y%m%d_%H%M%S")
    attempt = 0
    snapshot_dir = history_root / timestamp
    while snapshot_dir.exists():
        attempt += 1
        snapshot_dir = history_root / f"{timestamp}_{attempt:02d}"
    snapshot_dir.mkdir(parents=True, exist_ok=True)

    copied_any = False
    try:
        if target_yaml.exists():
            shutil.copy2(target_yaml, snapshot_dir / "target.yaml")
            copied_any = True
        if prep_dir.exists():
            shutil.copytree(prep_dir, snapshot_dir / "prep")
            copied_any = True
    except Exception:
        # Cleanup partially created snapshot on failure
        shutil.rmtree(snapshot_dir, ignore_errors=True)
        raise

    if not copied_any:
        shutil.rmtree(snapshot_dir, ignore_errors=True)
        return None

    _prune_old_snapshots(history_root, keep=keep)
    return snapshot_dir


def _prune_old_snapshots(history_root: Path, *, keep: int = 8) -> None:
    if keep <= 0 or not history_root.exists():
        return
    snapshots = [p for p in history_root.iterdir() if p.is_dir()]
    if len(snapshots) <= keep:
        return
    snapshots.sort(key=lambda p: p.stat().st_mtime)
    for candidate in snapshots[:-keep]:
        shutil.rmtree(candidate, ignore_errors=True)


__all__ = ["PipelineError", "run_manage_rfa", "init_decide_prep"]
