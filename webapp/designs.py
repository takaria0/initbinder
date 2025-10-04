"""Design pipeline orchestration helpers."""

from __future__ import annotations

import json
import math
import re
import shlex
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

import yaml

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
) -> List[str]:
    args: List[str] = [request.pdb_id]
    for arm in arms:
        args.extend(["--arm", arm])
    args.extend(["--total", str(request.total_designs)])
    args.extend(["--designs_per_task", str(designs_per_task)])
    args.extend(["--num_seq", str(request.num_sequences)])
    args.extend(["--temp", str(request.temperature)])
    if request.binder_chain_id:
        args.extend(["--binder_chain_id", request.binder_chain_id])
    if run_label:
        args.extend(["--run_tag", run_label])
    return args


def run_design_workflow(request: DesignRunRequest, *, job_store: JobStore, job_id: str) -> None:
    job_store.update(job_id, status=JobStatus.RUNNING, message="Generating pipeline scripts")
    log_buffer: List[str] = []

    cfg = load_config()
    workspace = cfg.paths.workspace_root or cfg.paths.project_root
    arms = _discover_arms(workspace, request.pdb_id)
    if not arms:
        raise ValueError("Failed to determine design arms; ensure prep-target has completed successfully.")

    run_label = request.run_label or datetime.utcnow().strftime("%Y%m%d_%H%M%S")
    designs_per_task = max(1, math.ceil(request.total_designs / max(len(arms), 1)))

    job_store.update(
        job_id,
        message=f"Generating pipeline scripts ({len(arms)} arms)",
        details={
            "arms": arms,
            "designs_per_task": designs_per_task,
            "run_label": run_label,
        },
    )
    job_store.append_log(job_id, f"[arms] {', '.join(arms)}")
    job_store.append_log(job_id, f"[designs_per_task] {designs_per_task}")

    def _log(line: str) -> None:
        log_buffer.append(line)
        job_store.append_log(job_id, line)

    args = _build_pipeline_args(request, arms, designs_per_task, run_label)
    run_manage_rfa("pipeline", args, log=_log)

    launcher = _parse_launcher_path(log_buffer)
    if not launcher or not launcher.exists():
        raise FileNotFoundError("Could not locate generated launcher script; check manage_rfa output")

    _instrument_launcher(launcher)
    job_store.update(job_id, message=f"Launcher instrumented: {launcher}")

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
        details={
            "job_ids": job_ids,
            "remote_launch": str(remote_path),
        },
    )


__all__ = ["run_design_workflow"]
