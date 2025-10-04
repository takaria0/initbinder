"""Design pipeline orchestration helpers."""

from __future__ import annotations

import json
import math
import re
import shlex
import textwrap
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
    job_store.update(job_id, status=JobStatus.RUNNING, message="Preparing pipeline configuration")

    cfg = load_config()
    workspace = cfg.paths.workspace_root or cfg.paths.project_root
    arms = _discover_arms(workspace, request.pdb_id)
    if not arms:
        raise ValueError("Failed to determine design arms; ensure prep-target has completed successfully.")

    run_label = request.run_label or datetime.utcnow().strftime("%Y%m%d_%H%M%S")
    designs_per_task = max(1, math.ceil(request.total_designs / max(len(arms), 1)))

    job_store.update(
        job_id,
        message=f"Syncing target to cluster ({len(arms)} arms)",
        details={
            "arms": arms,
            "designs_per_task": designs_per_task,
            "run_label": run_label,
        },
    )
    job_store.append_log(job_id, f"[arms] {', '.join(arms)}")
    job_store.append_log(job_id, f"[designs_per_task] {designs_per_task}")

    cluster = ClusterClient()
    sync_result, backup_rel = cluster.sync_target(request.pdb_id)
    if backup_rel:
        job_store.append_log(job_id, f"[sync] Remote target backup → {backup_rel}")
    if sync_result.stdout:
        for line in sync_result.stdout.splitlines():
            job_store.append_log(job_id, line)
    if sync_result.stderr:
        err = sync_result.stderr.strip()
        if err:
            job_store.append_log(job_id, err)

    job_store.update(job_id, message="Syncing tools to cluster")
    tools_result = cluster.sync_tools()
    if tools_result.stdout:
        for line in tools_result.stdout.splitlines():
            job_store.append_log(job_id, line)
    if tools_result.stderr:
        err = tools_result.stderr.strip()
        if err:
            job_store.append_log(job_id, err)

    pipeline_args = _build_pipeline_args(request, arms, designs_per_task, run_label)
    remote_cmd = shlex.join(["python", "manage_rfa.py", "pipeline", *pipeline_args])
    job_store.update(job_id, message="Generating pipeline scripts on cluster")
    pipeline_result = cluster.run(remote_cmd)

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

    instrument_script = textwrap.dedent(
        f"""
        python - <<'PY'
from pathlib import Path
import re

path = Path({launcher_path!r})
lines = path.read_text().splitlines()
assign_pattern = re.compile(r'^([A-Za-z0-9_]+)=\$\(\s*sbatch\\b')
patched = []
for line in lines:
    patched.append(line)
    m = assign_pattern.match(line.strip())
    if m:
        var_name = m.group(1)
        patched.append(f'echo "[JOB] {{var_name}}=${{{{var_name}}}}"')
path.write_text("\n".join(patched) + "\n")
PY
        """
    ).strip()
    cluster.run(instrument_script, check=True)
    job_store.append_log(job_id, f"[instrument] Patched launcher {launcher_path}")

    job_store.update(job_id, message="Submitting SLURM pipeline")
    sbatch_result = cluster.run(f"bash {shlex.quote(launcher_path)}")
    job_ids: Dict[str, str] = {}
    if sbatch_result.stdout:
        for line in sbatch_result.stdout.splitlines():
            job_store.append_log(job_id, line)
            match = _JOB_ECHO_RE.match(line.strip())
            if match:
                job_ids[match.group("var")] = match.group("jid")
    if sbatch_result.stderr:
        err = sbatch_result.stderr.strip()
        if err:
            job_store.append_log(job_id, err)

    stage2_ids = [jid for key, jid in job_ids.items() if key.startswith("jid_af3s2") and jid]
    job_store.update(
        job_id,
        details={
            "job_ids": job_ids,
            "remote_launch": launcher_path,
            "assessment_dependencies": stage2_ids,
        },
    )

    if stage2_ids:
        assessment_job_id = cluster.submit_assessment(
            pdb_id=request.pdb_id,
            binder_chain=request.binder_chain_id or "H",
            run_label=run_label,
            dependencies=stage2_ids,
            include_keyword=run_label,
        )
        if assessment_job_id:
            job_store.append_log(job_id, f"[assessment] Scheduled sbatch job {assessment_job_id}")
            job_store.update(job_id, details={"assessment_job_id": assessment_job_id})

    job_store.update(
        job_id,
        status=JobStatus.SUCCESS,
        message="Pipeline submitted to cluster",
    )


__all__ = ["run_design_workflow"]
