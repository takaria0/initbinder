"""Design pipeline orchestration helpers with pluggable engines."""

from __future__ import annotations

import json
import math
import os
import re
import shlex
import threading
import time
from abc import ABC, abstractmethod
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import yaml

from .config import load_config
from .hpc import ClusterClient
from .job_store import JobStatus, JobStore
from .models import DesignRunRequest


@dataclass(frozen=True)
class DesignEngineField:
    """UI metadata describing how a request field should be presented."""

    field_id: str
    label: str
    description: Optional[str] = None
    visible: bool = True
    debug_only: bool = False


@dataclass(frozen=True)
class DesignEngineMetadata:
    """Describes a concrete model engine exposed through the design workflow."""

    engine_id: str
    label: str
    description: str
    is_default: bool = False
    ui_fields: List[DesignEngineField] = field(default_factory=list)


class DesignEngine(ABC):
    """Abstract base for binder design engines."""

    metadata: DesignEngineMetadata

    @abstractmethod
    def run(
        self,
        request: DesignRunRequest,
        *,
        job_store: JobStore,
        job_id: str,
    ) -> None:
        """Execute the engine-specific workflow."""

    def to_metadata(self) -> DesignEngineMetadata:
        return self.metadata


class DesignEngineRegistry:
    """Registry and resolver for available design engines."""

    def __init__(self) -> None:
        self._engines: Dict[str, DesignEngine] = {}
        self._default_id: Optional[str] = None

    def register(self, engine: DesignEngine) -> None:
        meta = engine.to_metadata()
        key = meta.engine_id
        if not key:
            raise ValueError("Design engine id cannot be empty")
        self._engines[key] = engine
        if meta.is_default or self._default_id is None:
            self._default_id = key

    def get(self, engine_id: str) -> Optional[DesignEngine]:
        return self._engines.get(engine_id)

    def list_metadata(self) -> List[DesignEngineMetadata]:
        return [engine.to_metadata() for engine in self._engines.values()]

    def default_engine_id(self) -> str:
        if not self._default_id:
            raise RuntimeError("No design engines registered")
        return self._default_id


_ENGINE_REGISTRY = DesignEngineRegistry()


@dataclass
class BoltzGenSpecInfo:
    arm: str
    spec_path: Path
    rel_spec_path: Path
    output_rel: Path
    num_designs: int
    hotspot_count: int


class RFAntibodyEngine(DesignEngine):
    """Reproduces the existing RFdiffusion → ProteinMPNN → AF3 workflow."""

    metadata = DesignEngineMetadata(
        engine_id="rfantibody",
        label="RFantibody (RFdiffusion → MPNN → AF3)",
        description="Legacy antibody pipeline using RFdiffusion, ProteinMPNN, and AlphaFold3.",
        is_default=True,
        ui_fields=[
            DesignEngineField(
                field_id="total_designs",
                label="Total designs",
                description="Split evenly across detected epitope arms.",
            ),
            DesignEngineField(
                field_id="num_sequences",
                label="Sequences per backbone",
                description="ProteinMPNN sequences to sample per RFdiffusion backbone.",
            ),
            DesignEngineField(
                field_id="temperature",
                label="RFdiffusion temperature",
                description="Controls sampling diversity for RFdiffusion (0.0–1.0).",
            ),
            DesignEngineField(
                field_id="binder_chain_id",
                label="Binder chain ID",
                description="Override default binder chain for downstream assessment (debug only).",
                debug_only=True,
            ),
            DesignEngineField(
                field_id="af3_seed",
                label="AlphaFold 3 seed",
                description="Seed forwarded to AlphaFold 3 stage for reproducibility.",
            ),
            DesignEngineField(
                field_id="rfdiff_crop_radius",
                label="RFdiffusion crop radius",
                description="Enable to crop the prepared target around hotspots for RFdiffusion.",
                debug_only=True,
            ),
        ],
    )

    _LAUNCH_PATH_RE = re.compile(r"\[ok] Wrote launcher: (?P<path>.*)")
    _JOB_ECHO_RE = re.compile(r"^\[JOB] (?P<var>[A-Za-z0-9_]+)=(?P<jid>\d+)")
    _SCRIPT_PATH_RE = re.compile(r"sbatch[^|;]*?([^\s]+\.sh)")

    def run(
        self,
        request: DesignRunRequest,
        *,
        job_store: JobStore,
        job_id: str,
    ) -> None:
        job_store.update(job_id, status=JobStatus.RUNNING, message="Preparing pipeline configuration")
        job_store.append_log(job_id, f"[engine] {self.metadata.engine_id}")

        cfg = load_config()
        workspace = cfg.paths.workspace_root or cfg.paths.project_root
        arms = self._discover_arms(workspace, request.pdb_id)
        if not arms:
            raise ValueError(
                "Failed to determine design arms; ensure prep-target has completed successfully."
            )

        run_label = request.run_label or datetime.utcnow().strftime("%Y%m%d_%H%M%S")
        designs_per_task = max(1, math.ceil(request.total_designs / max(len(arms), 1)))
        binder_chain = (request.binder_chain_id or "H").strip().upper() or "H"
        assessment_shards = max(1, math.ceil(request.total_designs / 1000))

        job_store.update(
            job_id,
            message=f"Syncing target to cluster ({len(arms)} arms)",
            details={
                "model_engine": self.metadata.engine_id,
                "arms": arms,
                "designs_per_task": designs_per_task,
                "run_label": run_label,
                "af3_seed": request.af3_seed,
                "binder_chain": binder_chain,
                "run_assess": request.run_assess,
                "assessment_shards": assessment_shards,
                "rfdiff_crop_radius": request.rfdiff_crop_radius,
            },
        )
        job_store.append_log(job_id, f"[arms] {', '.join(arms)}")
        job_store.append_log(job_id, f"[designs_per_task] {designs_per_task}")
        job_store.append_log(job_id, f"[af3_seed] {request.af3_seed}")
        job_store.append_log(job_id, f"[binder_chain] {binder_chain}")
        job_store.append_log(job_id, f"[run_assess] {request.run_assess}")
        if request.rfdiff_crop_radius is None or request.rfdiff_crop_radius <= 0:
            job_store.append_log(job_id, "[rfdiff_crop] disabled")
        else:
            job_store.append_log(
                job_id,
                f"[rfdiff_crop] radius={request.rfdiff_crop_radius}Å",
            )
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

        pipeline_args = self._build_pipeline_args(request, arms, designs_per_task, run_label, binder_chain)
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
        remote_cmd = (f"{env_prefix} " if env_prefix else "") + shlex.join(
            [
                "python",
                "manage_rfa.py",
                "pipeline",
                *pipeline_args,
            ]
        )
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
            launcher = self._parse_launcher_path(log_buffer)
            if launcher:
                launcher_path = str(launcher)
        if not launcher_path:
            raise FileNotFoundError(
                "Could not locate generated launcher script on the cluster; check manage_rfa output"
            )

        cat_cmd = f"cat {shlex.quote(launcher_path)}"
        launcher_contents = cluster.run(cat_cmd).stdout or ""
        script_order: List[str] = []
        for line in launcher_contents.splitlines():
            match = self._SCRIPT_PATH_RE.search(line)
            if match:
                script_order.append(match.group(1))

        job_store.update(job_id, message="Submitting SLURM pipeline")
        submit_cmd = (
            "conda deactivate >/dev/null 2>&1 || true; "
            "conda deactivate >/dev/null 2>&1 || true; "
            f"bash {shlex.quote(launcher_path)}"
        )
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
                "model_engine": self.metadata.engine_id,
            },
        )

        tracked_for_monitor: List[str] = list(stage2_ids)
        if assessment_job_id:
            tracked_for_monitor.append(assessment_job_id)
        if assessment_merge_job_id:
            tracked_for_monitor.append(assessment_merge_job_id)
        self._start_squeue_monitor(cluster, job_store, job_id, tracked_for_monitor)

        if request.run_assess and pipeline_assess_job and stage2_ids:
            # If all stage-2 jobs stall on dependency reasons the rescue monitor
            # will submit a fresh assessment job with no dependencies so users can
            # still review partial results.
            self._schedule_assessment_rescue(
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

    @staticmethod
    def _parse_launcher_path(log_lines: List[str]) -> Optional[Path]:
        for line in reversed(log_lines):
            match = RFAntibodyEngine._LAUNCH_PATH_RE.search(line)
            if match:
                return Path(match.group("path")).expanduser().resolve()
        return None

    @staticmethod
    def _sanitize_epitope_name(name: str) -> str:
        s = str(name).strip()
        return s.replace(" ", "_").replace("/", "_").replace("\\", "_")

    @staticmethod
    def _load_epitope_names(workspace: Path, pdb_id: str) -> List[str]:
        prep_dir = workspace / "targets" / pdb_id.upper() / "prep"
        metadata_path = prep_dir / "epitopes_metadata.json"
        names: List[str] = []
        if metadata_path.exists():
            try:
                data = json.loads(metadata_path.read_text()) or {}
                names = [
                    str(ep.get("name")).strip()
                    for ep in data.get("epitopes", [])
                    if ep.get("name")
                ]
            except Exception:
                names = []
        if not names:
            target_yaml = workspace / "targets" / pdb_id.upper() / "target.yaml"
            if target_yaml.exists():
                try:
                    cfg = yaml.safe_load(target_yaml.read_text()) or {}
                    names = [
                        str(ep.get("name")).strip()
                        for ep in (cfg.get("epitopes") or [])
                        if ep.get("name")
                    ]
                except Exception:
                    names = []
        return [name for name in names if name]

    @staticmethod
    def _detect_hotspot_variants(prep_dir: Path, epitope_name: str) -> List[str]:
        san = RFAntibodyEngine._sanitize_epitope_name(epitope_name)
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

    @classmethod
    def _discover_arms(cls, workspace: Path, pdb_id: str) -> List[str]:
        prep_dir = workspace / "targets" / pdb_id.upper() / "prep"
        if not prep_dir.exists():
            raise FileNotFoundError(
                f"Prep directory not found for {pdb_id.upper()} (expected {prep_dir})"
            )

        epitope_names = cls._load_epitope_names(workspace, pdb_id)
        if not epitope_names:
            raise ValueError("No epitopes available; run prep-target before submitting designs.")

        arms: List[str] = []
        for name in epitope_names:
            variants = cls._detect_hotspot_variants(prep_dir, name)
            for variant in variants:
                arms.append(f"{name}@{variant}")
        return list(dict.fromkeys(arms))

    @staticmethod
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
        if request.rfdiff_crop_radius is not None:
            args.extend(["--crop_radius", str(request.rfdiff_crop_radius)])
        return args

    @staticmethod
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
            or "user"
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

    @classmethod
    def _schedule_assessment_rescue(
        cls,
        cluster: ClusterClient,
        job_store: JobStore,
        job_id: str,
        *,
        stage2_ids: List[str],
        pdb_id: str,
        binder_chain: str,
        run_label: str,
    ) -> None:
        """Submit an assessment job without dependencies if stage2 jobs are stuck."""

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
                    cls._start_squeue_monitor(cluster, job_store, job_id, [rescue_id])
                return

        threading.Thread(
            target=_poll,
            name=f"assess-rescue-{job_id}",
            daemon=True,
        ).start()

class BoltzGenEngine(DesignEngine):
    """BoltzGen diffusion-based binder design pipeline."""

    metadata = DesignEngineMetadata(
        engine_id="boltzgen",
        label="BoltzGen Diffusion",
        description="BoltzGen generative diffusion pipeline for protein binder design.",
        ui_fields=[
            DesignEngineField(
                field_id="total_designs",
                label="Total designs",
                description="Total BoltzGen designs (single spec; no epitope splitting).",
            ),
            DesignEngineField(
                field_id="num_sequences",
                label="Sequences per backbone",
                visible=False,
            ),
            DesignEngineField(
                field_id="temperature",
                label="RFdiffusion temperature",
                visible=False,
            ),
            DesignEngineField(
                field_id="binder_chain_id",
                label="Binder chain ID",
                visible=False,
                debug_only=True,
            ),
            DesignEngineField(
                field_id="af3_seed",
                label="AlphaFold 3 seed",
                visible=False,
            ),
            DesignEngineField(
                field_id="rfdiff_crop_radius",
                label="RFdiffusion crop radius",
                visible=False,
                debug_only=True,
            ),
            DesignEngineField(
                field_id="boltzgen_crop_radius",
                label="BoltzGen crop radius",
                description="Crop target around hotspots (Å) for BoltzGen configs/specs.",
            ),
        ],
    )

    DEFAULT_SEQUENCE_MIN = 85
    DEFAULT_SEQUENCE_MAX = 115
    DEFAULT_PARTITION = "gpu"
    DEFAULT_GPUS = "A100:1"
    DEFAULT_CPUS = 8
    DEFAULT_MEM_GB = 64
    DEFAULT_TIME_H = 12

    _HOTSPOT_RE = re.compile(r"^([A-Za-z0-9]+)(\d+)$")
    _RES_TOKEN_DELIM_RE = re.compile(
        r"^\s*(?P<chain>[A-Za-z0-9]+)\s*[:_\-]\s*(?P<resnum>-?\d+)\s*(?P<icode>[A-Za-z]?)\s*$"
    )
    _RES_TOKEN_PLAIN_RE = re.compile(
        r"^\s*(?P<chain>[A-Za-z0-9]+?)\s*(?P<resnum>-?\d+)\s*(?P<icode>[A-Za-z]?)\s*$"
    )

    def run(
        self,
        request: DesignRunRequest,
        *,
        job_store: JobStore,
        job_id: str,
    ) -> None:
        job_store.update(job_id, status=JobStatus.RUNNING, message="Preparing BoltzGen specifications")
        job_store.append_log(job_id, f"[engine] {self.metadata.engine_id}")

        cfg = load_config()
        workspace = cfg.paths.workspace_root or cfg.paths.project_root
        run_label = request.run_label or datetime.utcnow().strftime("%Y%m%d_%H%M%S")
        try:
            arms = RFAntibodyEngine._discover_arms(workspace, request.pdb_id)
        except Exception as exc:
            arms = []
            job_store.append_log(
                job_id,
                f"[boltzgen] warn: unable to enumerate prep arms ({exc}); proceeding without epitope constraints.",
            )
        total_designs = max(1, request.total_designs)
        if arms:
            job_store.append_log(
                job_id,
                f"[boltzgen] detected {len(arms)} prep arm(s); generating a single full-target spec without epitope filters.",
            )
        else:
            job_store.append_log(
                job_id,
                "[boltzgen] no prep arms detected; generating a full-target BoltzGen spec.",
            )

        spec_infos = self._prepare_specs(
            workspace=workspace,
            pdb_id=request.pdb_id,
            run_label=run_label,
            num_designs=total_designs,
            binding_override=request.boltz_binding,
            crop_radius=request.boltzgen_crop_radius,
        )
        details = {
            "model_engine": self.metadata.engine_id,
            "arms": arms,
            "run_label": run_label,
            "num_designs": total_designs,
            "spec_paths": [str(info.rel_spec_path) for info in spec_infos],
            "spec_outputs": [str(info.output_rel) for info in spec_infos],
        }
        job_store.update(job_id, message=f"Syncing target to cluster ({len(spec_infos)} specs)", details=details)
        for info in spec_infos:
            job_store.append_log(
                job_id,
                f"[spec] {info.arm} -> {info.rel_spec_path} (hotspots={info.hotspot_count})",
            )
        job_store.append_log(
            job_id,
            f"[boltzgen] binder sequence length {self.DEFAULT_SEQUENCE_MIN}..{self.DEFAULT_SEQUENCE_MAX} residues",
        )

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

        binder_root = cluster.remote_root
        if binder_root is None:
            raise RuntimeError("Cluster remote_root not configured; cannot launch BoltzGen pipeline")
        remote_root_path = Path(binder_root)

        env_vars = {"INITBINDER_ROOT": remote_root_path}
        target_base: Optional[Path] = None
        if cluster.target_root:
            target_base = Path(cluster.target_root)
        elif cluster.remote_root:
            target_base = Path(cluster.remote_root)
        if target_base is not None:
            if target_base.name.lower() != "targets":
                target_base = target_base / "targets"
            env_vars["INITBINDER_TARGET_ROOT"] = target_base

        bg_cfg = cluster.cfg.boltzgen
        partition = bg_cfg.partition or self.DEFAULT_PARTITION
        gpus = bg_cfg.gpus or self.DEFAULT_GPUS
        cpus = int(bg_cfg.cpus or self.DEFAULT_CPUS)
        mem_gb = int(bg_cfg.mem_gb or self.DEFAULT_MEM_GB)
        time_h = int(request.boltz_time_hours or bg_cfg.time_hours or self.DEFAULT_TIME_H)
        job_store.append_log(
            job_id,
            "[boltzgen] resources partition=%s account=%s gpus=%s cpus=%s mem=%sG time=%sh"
            % (partition, bg_cfg.account or "-", gpus, cpus, mem_gb, time_h),
        )
        if bg_cfg.output_root:
            output_root_path = Path(bg_cfg.output_root)
        else:
            output_root_path = Path("targets") / request.pdb_id.upper() / "designs" / "_boltzgen" / run_label
        remote_output_root = self._resolve_remote_workspace_path(
            output_root_path, remote_root=remote_root_path, target_base=target_base
        )
        output_root_value = str(remote_output_root)
        remote_cmd_parts: List[str] = [
            "python",
            "tools/boltzgen/pipeline.py",
            "pipeline",
            request.pdb_id,
            "--run_label",
            run_label,
            "--num_designs",
            str(total_designs),
            "--protocol",
            self._resolve_protocol(bg_cfg),
            "--scripts_dir",
            "tools/boltzgen",
            "--launcher_dir",
            "tools/launchers",
            "--output_root",
            output_root_value,
            "--partition",
            partition,
            "--gpus",
            gpus,
            "--cpus",
            str(cpus),
            "--mem",
            f"{mem_gb}G",
            "--time_h",
            str(time_h),
        ]
        if bg_cfg.account:
            remote_cmd_parts.extend(["--account", bg_cfg.account])
        conda_cmd = bg_cfg.conda_activate or cluster.cfg.conda_activate
        if conda_cmd:
            remote_cmd_parts.extend(["--conda_activate", conda_cmd])
        if bg_cfg.cache_dir:
            remote_cmd_parts.extend(["--cache_dir", str(bg_cfg.cache_dir)])
        for info in spec_infos:
            remote_spec_path = self._resolve_remote_workspace_path(
                info.rel_spec_path, remote_root=remote_root_path, target_base=target_base
            )
            remote_cmd_parts.extend(["--spec", str(remote_spec_path)])
        remote_cmd_parts.append("--submit")
        extra_args = list(bg_cfg.extra_run_args or [])
        if extra_args:
            remote_cmd_parts.append("--extra_run_args")
            remote_cmd_parts.extend(extra_args)

        env_prefix = " ".join(f"{key}={shlex.quote(str(value))}" for key, value in env_vars.items())
        remote_cmd = (f"{env_prefix} " if env_prefix else "") + shlex.join(remote_cmd_parts)
        job_store.update(job_id, message="Generating BoltzGen pipeline scripts on cluster")
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

        launcher_path, job_ids, returned_label = self._parse_pipeline_output(log_buffer)
        if launcher_path is None:
            raise FileNotFoundError(
                "BoltzGen pipeline did not report a launcher script; ensure tools/boltzgen/pipeline.py ran successfully."
            )
        if returned_label and returned_label != run_label:
            job_store.append_log(
                job_id,
                f"[warn] Pipeline reported run label {returned_label}; expected {run_label}",
            )

        job_store.update(
            job_id,
            details={
                **details,
                "remote_launch": launcher_path,
                "job_ids": job_ids,
                "boltzgen_conda": conda_cmd,
                "boltzgen_partition": partition,
                "boltzgen_account": bg_cfg.account,
                "boltzgen_gpus": gpus,
                "boltzgen_time_h": time_h,
                "run_assess_supported": False,
                "boltzgen_output_root": output_root_value,
                "boltzgen_cache_dir": str(bg_cfg.cache_dir) if bg_cfg.cache_dir else None,
            },
        )
        if job_ids:
            RFAntibodyEngine._start_squeue_monitor(cluster, job_store, job_id, job_ids)

        if request.run_assess:
            job_store.append_log(
                job_id,
                "[warn] Automatic assessment is not yet supported for BoltzGen runs; skipping run_assess.",
            )

        job_store.update(
            job_id,
            status=JobStatus.SUCCESS,
            message="BoltzGen pipeline submitted to cluster",
        )

    @staticmethod
    def _load_epitope_metadata(prep_dir: Path) -> List[dict]:
        meta_path = prep_dir / "epitopes_metadata.json"
        if not meta_path.exists():
            return []
        try:
            data = json.loads(meta_path.read_text())
            epitopes = data.get("epitopes") or []
            output: List[dict] = []
            for ep in epitopes:
                name = (ep.get("name") or "").strip()
                if not name:
                    continue
                mask = ep.get("files", {}).get("mask_json")
                mask_residues: List[str] = []
                if mask:
                    mask_path = prep_dir / mask
                    if mask_path.exists():
                        try:
                            mask_data = json.loads(mask_path.read_text())
                            if isinstance(mask_data, list):
                                mask_residues = [str(x).strip() for x in mask_data if str(x).strip()]
                        except Exception:
                            mask_residues = []
                hotspots = ep.get("hotspots") or []
                hotspots = [str(h).strip() for h in hotspots if str(h).strip()]
                output.append(
                    {
                        "name": name,
                        "hotspots": hotspots,
                        "mask_residues": mask_residues,
                    }
                )
            return output
        except Exception:
            return []

    @staticmethod
    def _distribute_designs(total: int, count: int) -> List[int]:
        if count <= 0:
            return []
        base = total // count
        remainder = total % count
        designs = [base + (1 if i < remainder else 0) for i in range(count)]
        # Ensure no zero-design specs; allow slight over-allocation if total < count.
        designs = [d if d > 0 else 1 for d in designs]
        return designs

    def _prepare_specs(
        self,
        *,
        workspace: Path,
        pdb_id: str,
        run_label: str,
        num_designs: int,
        binding_override: Optional[str] = None,
        crop_radius: Optional[float] = None,
    ) -> List[BoltzGenSpecInfo]:
        target_dir = workspace / "targets" / pdb_id.upper()
        prep_dir = target_dir / "prep"
        if not prep_dir.exists():
            raise FileNotFoundError(f"Prep directory not found for {pdb_id.upper()} (expected {prep_dir})")
        prepared_pdb = prep_dir / "prepared.pdb"
        if not prepared_pdb.exists():
            raise FileNotFoundError(f"Missing prepared.pdb for {pdb_id.upper()} (expected {prepared_pdb})")

        spec_root = target_dir / "designs" / "_boltzgen" / run_label
        spec_root.mkdir(parents=True, exist_ok=True)

        spec_filename = "full_target.yaml"
        spec_path = spec_root / spec_filename
        scaffold_paths = self._resolve_nanobody_scaffolds()

        binding_override_map: Optional[Dict[str, List[int]]] = None
        if binding_override:
            binding_override_map = self._expand_epitope_residues([binding_override])
        ep_meta = self._load_epitope_metadata(prep_dir) if not binding_override_map else []

        if ep_meta:
            specs: List[BoltzGenSpecInfo] = []
            designs_split = self._distribute_designs(num_designs, len(ep_meta))
            for idx, ep in enumerate(ep_meta):
                san = RFAntibodyEngine._sanitize_epitope_name(ep["name"])
                ep_spec = spec_root / f"{san or f'epitope{idx+1}'}.yaml"
                info = self._write_spec(
                    workspace=workspace,
                    spec_path=ep_spec,
                    prepared_pdb=prepared_pdb,
                    pdb_id=pdb_id,
                    run_label=run_label,
                    arm=ep["name"],
                    hotspot_keys=ep.get("hotspots") or [],
                    epitope_residues=ep.get("mask_residues") or [],
                    num_designs=designs_split[idx],
                    scaffold_paths=scaffold_paths,
                    binding_override=None,
                    crop_radius=crop_radius,
                )
                specs.append(info)
            return specs

        info = self._write_spec(
            workspace=workspace,
            spec_path=spec_path,
            prepared_pdb=prepared_pdb,
            pdb_id=pdb_id,
            run_label=run_label,
            arm="FULL_TARGET",
            hotspot_keys=[],
            epitope_residues=[],
            num_designs=num_designs,
            scaffold_paths=scaffold_paths,
            binding_override=binding_override_map,
            crop_radius=crop_radius,
        )
        return [info]

    @staticmethod
    def _resolve_protocol(bg_cfg) -> str:
        value = getattr(bg_cfg, "protocol", None)
        if value:
            value = str(value).strip()
        return value or "protein-anything"

    @staticmethod
    def _resolve_nanobody_scaffolds() -> List[str]:
        cfg = load_config()
        bg_cfg = getattr(cfg.cluster, "boltzgen", None)
        if not bg_cfg:
            return []
        candidates = getattr(bg_cfg, "nanobody_scaffolds", None) or []
        if not isinstance(candidates, list):
            return []
        return [str(entry).strip() for entry in candidates if str(entry).strip()]

    @staticmethod
    def _resolve_remote_workspace_path(
        path: Path | str,
        *,
        remote_root: Path,
        target_base: Optional[Path],
    ) -> Path:
        candidate = Path(path)
        if candidate.is_absolute():
            return candidate
        parts = candidate.parts
        if target_base and parts and parts[0].lower() == "targets":
            remainder = Path(*parts[1:]) if len(parts) > 1 else Path()
            return target_base / remainder if remainder.parts else target_base
        return remote_root / candidate

    @staticmethod
    def _load_epitope_residue_map(target_yaml: Path) -> Dict[str, List[str]]:
        if not target_yaml.exists():
            return {}
        try:
            data = yaml.safe_load(target_yaml.read_text()) or {}
        except Exception:
            return {}
        epitope_map: Dict[str, List[str]] = {}
        for entry in data.get("epitopes", []) or []:
            name = str(entry.get("name") or "").strip()
            residues = [str(r).strip() for r in entry.get("residues", []) or [] if str(r).strip()]
            if name:
                epitope_map[name] = residues
        return epitope_map

    @staticmethod
    def _expand_epitope_residues(tokens: Iterable[str]) -> Dict[str, List[int]]:
        mapping: Dict[str, set[int]] = defaultdict(set)
        for token in tokens:
            token = (token or "").strip()
            if not token or ":" not in token:
                continue
            chain, rest = token.split(":", 1)
            parts = [part.strip() for part in rest.replace("..", "-").split(",") if part.strip()]
            for part in parts:
                if "-" in part:
                    start_str, end_str = part.split("-", 1)
                    start = int(start_str)
                    end = int(end_str)
                    lo, hi = sorted((start, end))
                    for res in range(lo, hi + 1):
                        mapping[chain].add(res)
                else:
                    mapping[chain].add(int(part))
        return {chain: sorted(values) for chain, values in mapping.items()}

    def _parse_hotspot_map(self, keys: Iterable[str]) -> Dict[str, List[int]]:
        mapping: Dict[str, List[int]] = defaultdict(list)
        for key in keys:
            parsed = self._parse_chain_res_token(key)
            if not parsed:
                continue
            chain, res = parsed
            mapping[chain].append(res)
        return {chain: sorted(set(vals)) for chain, vals in mapping.items()}

    @staticmethod
    def _parse_chain_res_token(token: object) -> Optional[Tuple[str, int]]:
        if token is None:
            return None
        text = str(token).strip()
        if not text:
            return None
        match = BoltzGenEngine._RES_TOKEN_DELIM_RE.match(text)
        if not match:
            match = BoltzGenEngine._RES_TOKEN_PLAIN_RE.match(text)
        if not match:
            return None
        chain = (match.group("chain") or "").strip()
        res_s = match.group("resnum")
        try:
            resnum = int(res_s)
        except Exception:
            return None
        if not chain:
            return None
        return chain, resnum

    @staticmethod
    def _collect_mmcif_residue_atoms(
        mmcif_path: Path,
        allowed_chains: Optional[set[str]] = None,
    ) -> Dict[Tuple[str, int], List[Tuple[float, float, float]]]:
        try:
            from Bio.PDB.MMCIF2Dict import MMCIF2Dict
        except ImportError:
            return {}
        try:
            mmcif = MMCIF2Dict(str(mmcif_path))
        except Exception:
            return {}
        chains = mmcif.get("_atom_site.label_asym_id") or []
        seqs = mmcif.get("_atom_site.label_seq_id") or []
        xs = mmcif.get("_atom_site.Cartn_x") or []
        ys = mmcif.get("_atom_site.Cartn_y") or []
        zs = mmcif.get("_atom_site.Cartn_z") or []
        elems = mmcif.get("_atom_site.type_symbol") or []
        if isinstance(chains, str):
            chains = [chains]
        if isinstance(seqs, str):
            seqs = [seqs]
        if isinstance(xs, str):
            xs = [xs]
        if isinstance(ys, str):
            ys = [ys]
        if isinstance(zs, str):
            zs = [zs]
        if isinstance(elems, str):
            elems = [elems]
        if not (chains and seqs and xs and ys and zs):
            return {}
        if not elems or len(elems) != len(chains):
            elems = [None] * len(chains)
        allowed_upper = {c.upper() for c in allowed_chains} if allowed_chains else None
        out: Dict[Tuple[str, int], List[Tuple[float, float, float]]] = defaultdict(list)
        for chain, seq, x, y, z, elem in zip(chains, seqs, xs, ys, zs, elems):
            chain_id = str(chain).strip()
            if not chain_id:
                continue
            if allowed_upper and chain_id.upper() not in allowed_upper:
                continue
            seq_txt = str(seq).strip()
            if seq_txt in {".", "?", ""}:
                continue
            try:
                seq_id = int(float(seq_txt))
            except Exception:
                continue
            elem_txt = (str(elem).strip().upper() if elem is not None else "")
            if elem_txt.startswith("H") or elem_txt == "D":
                continue
            try:
                coord = (float(x), float(y), float(z))
            except Exception:
                continue
            out[(chain_id, seq_id)].append(coord)
        return out

    @staticmethod
    def _collect_pdb_residue_atoms(
        pdb_path: Path,
        allowed_chains: Optional[set[str]] = None,
    ) -> Dict[Tuple[str, int], List[Tuple[float, float, float]]]:
        try:
            from Bio.PDB import PDBParser
        except ImportError:
            return {}
        try:
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure("boltz", str(pdb_path))
        except Exception:
            return {}
        allowed_upper = {c.upper() for c in allowed_chains} if allowed_chains else None
        out: Dict[Tuple[str, int], List[Tuple[float, float, float]]] = defaultdict(list)
        try:
            model = structure[0]
        except Exception:
            return {}
        for chain in model:
            chain_id = str(chain.id).strip()
            if not chain_id:
                continue
            if allowed_upper and chain_id.upper() not in allowed_upper:
                continue
            for residue in chain:
                if residue.id[0] != " ":
                    continue
                resnum = residue.id[1]
                coords: List[Tuple[float, float, float]] = []
                for atom in residue.get_atoms():
                    elem = str(getattr(atom, "element", "")).strip().upper()
                    if elem.startswith("H") or elem == "D":
                        continue
                    try:
                        coord = atom.coord
                        coords.append((float(coord[0]), float(coord[1]), float(coord[2])))
                    except Exception:
                        continue
                if coords:
                    out[(chain_id, int(resnum))].extend(coords)
        return out

    def _crop_res_index_by_hotspots(
        self,
        prepared_pdb: Path,
        hotspot_keys: Iterable[str],
        radius: float,
        allowed_chains: Optional[Sequence[str]] = None,
    ) -> Dict[str, List[int]]:
        if radius <= 0:
            return {}
        allowed = {str(c).strip() for c in (allowed_chains or []) if str(c).strip()}
        allowed = allowed or None
        if not prepared_pdb.exists():
            return {}
        suffix = prepared_pdb.suffix.lower()
        if suffix in {".cif", ".mmcif"}:
            residue_atoms = self._collect_mmcif_residue_atoms(prepared_pdb, allowed_chains=allowed)
        else:
            residue_atoms = self._collect_pdb_residue_atoms(prepared_pdb, allowed_chains=allowed)
        if not residue_atoms:
            return {}

        hotspot_coords: List[Tuple[float, float, float]] = []
        for key in hotspot_keys or []:
            parsed = self._parse_chain_res_token(key)
            if not parsed:
                continue
            chain, resnum = parsed
            coords = residue_atoms.get((chain, resnum))
            if coords:
                hotspot_coords.extend(coords)
        if not hotspot_coords:
            return {}

        radius2 = float(radius) ** 2
        selected: set[Tuple[str, int]] = set()
        for resid, coords in residue_atoms.items():
            if resid in selected:
                continue
            for coord in coords:
                cx, cy, cz = coord
                for hx, hy, hz in hotspot_coords:
                    dx = cx - hx
                    dy = cy - hy
                    dz = cz - hz
                    if (dx * dx + dy * dy + dz * dz) <= radius2:
                        selected.add(resid)
                        break
                if resid in selected:
                    break

        out: Dict[str, List[int]] = defaultdict(list)
        for chain, resnum in selected:
            out[chain].append(int(resnum))
        return {chain: sorted(set(resnums)) for chain, resnums in out.items()}

    @staticmethod
    def _format_ranges(numbers: Iterable[int]) -> str:
        seq = sorted(set(int(n) for n in numbers))
        if not seq:
            return ""
        ranges: List[tuple[int, int]] = []
        start = prev = seq[0]
        for n in seq[1:]:
            if n == prev + 1:
                prev = n
                continue
            ranges.append((start, prev))
            start = prev = n
        ranges.append((start, prev))
        parts = []
        for lo, hi in ranges:
            parts.append(str(lo) if lo == hi else f"{lo}..{hi}")
        return ",".join(parts)

    def _write_spec(
        self,
        *,
        workspace: Path,
        spec_path: Path,
        prepared_pdb: Path,
        pdb_id: str,
        run_label: str,
        arm: str,
        hotspot_keys: List[str],
        epitope_residues: List[str],
        num_designs: int,
        scaffold_paths: Optional[List[str]] = None,
        binding_override: Optional[Dict[str, List[int]]] = None,
        crop_radius: Optional[float] = None,
        crop_chains: Optional[Sequence[str]] = None,
    ) -> BoltzGenSpecInfo:
        spec_path.parent.mkdir(parents=True, exist_ok=True)
        hotspot_map = self._parse_hotspot_map(hotspot_keys)
        include_map = self._expand_epitope_residues(epitope_residues)
        use_res_index = False
        if crop_radius is not None and crop_radius > 0 and hotspot_keys:
            crop_map = self._crop_res_index_by_hotspots(
                prepared_pdb,
                hotspot_keys,
                crop_radius,
                allowed_chains=crop_chains,
            )
            if crop_map:
                include_map = crop_map
                use_res_index = True
        binding_map = binding_override if binding_override else (hotspot_map if hotspot_map else include_map)
        if binding_override and not include_map:
            include_map = binding_override

        rel_prepared = os.path.relpath(prepared_pdb, spec_path.parent)

        include_entries: List[Dict[str, Dict[str, str]]] = []
        for chain, residues in sorted(include_map.items()):
            entry = {"chain": {"id": chain}}
            if use_res_index:
                formatted = self._format_ranges(residues)
                if formatted:
                    entry["chain"]["res_index"] = formatted
            include_entries.append(entry)

        binding_entries: List[Dict[str, Dict[str, str]]] = []
        for chain, residues in sorted(binding_map.items()):
            formatted = self._format_ranges(residues)
            if formatted:
                binding_entries.append({"chain": {"id": chain, "binding": formatted}})

        file_block: Dict[str, object] = {"path": rel_prepared}
        if include_entries:
            file_block["include"] = include_entries
        if binding_entries:
            file_block["binding_types"] = binding_entries

        entities: List[Dict[str, object]] = [{"file": file_block}]
        scaffold_entries = [str(path).strip() for path in (scaffold_paths or []) if str(path).strip()]
        if scaffold_entries:
            entities.append({"file": {"path": scaffold_entries}})
        else:
            entities.append(
                {
                    "protein": {
                        "id": "DES",
                        "sequence": f"{self.DEFAULT_SEQUENCE_MIN}..{self.DEFAULT_SEQUENCE_MAX}",
                    }
                }
            )

        spec_data: Dict[str, object] = {"entities": entities}

        spec_path.write_text(yaml.safe_dump(spec_data, sort_keys=False))
        rel_spec = spec_path.relative_to(workspace)
        output_rel = (
            Path("targets")
            / pdb_id.upper()
            / "designs"
            / "_boltzgen"
            / run_label
            / spec_path.stem
        )
        hotspot_count = sum(len(res) for res in binding_map.values())
        return BoltzGenSpecInfo(
            arm=arm,
            spec_path=spec_path,
            rel_spec_path=rel_spec,
            output_rel=output_rel,
            num_designs=num_designs,
            hotspot_count=hotspot_count,
        )

    @staticmethod
    def _parse_pipeline_output(log_lines: List[str]) -> tuple[Optional[str], List[str], Optional[str]]:
        launcher: Optional[str] = None
        job_ids: List[str] = []
        run_label: Optional[str] = None
        for line in log_lines:
            if "[ok] Wrote launcher:" in line:
                launcher = line.split(":", 1)[1].strip()
            elif "Submitted batch job" in line:
                parts = line.strip().split()
                if parts:
                    job_ids.append(parts[-1])
            elif "[ok] Run label:" in line:
                run_label = line.split(":", 1)[1].strip()
        return launcher, job_ids, run_label


_ENGINE_REGISTRY.register(RFAntibodyEngine())
_ENGINE_REGISTRY.register(BoltzGenEngine())


def list_design_engines() -> List[DesignEngineMetadata]:
    """Return metadata for every registered design engine."""

    return _ENGINE_REGISTRY.list_metadata()


def get_design_engine(engine_id: str) -> Optional[DesignEngine]:
    """Resolve a design engine by id."""

    return _ENGINE_REGISTRY.get(engine_id)


def default_design_engine_id() -> str:
    """Return the registry's default engine id."""

    return _ENGINE_REGISTRY.default_engine_id()


def run_design_workflow(request: DesignRunRequest, *, job_store: JobStore, job_id: str) -> None:
    """Entry point invoked by the web workflows module."""

    engine_id = (request.model_engine or "").strip() or default_design_engine_id()
    engine = get_design_engine(engine_id)
    if engine is None:
        raise ValueError(f"Unknown design engine '{request.model_engine}'")
    engine.run(request, job_store=job_store, job_id=job_id)


__all__ = [
    "run_design_workflow",
    "list_design_engines",
    "get_design_engine",
    "default_design_engine_id",
    "DesignEngineMetadata",
    "DesignEngine",
]
